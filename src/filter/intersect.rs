use std::{cell::RefCell, cmp::{Ordering, Reverse}, collections::BinaryHeap, f64::consts::E, marker::PhantomData, mem, rc::Rc, slice};

use exact_number::BasedExpr;
use nalgebra::{DMatrix, RawStorage, RealField, Vector2};
use num_traits::{RefNum, Zero};

use crate::{filter::Coordinate, geom::{self, FloatOrd, LineIntersection, NumEx, SegmentIntersection, VectorView2Dyn}, Frame};

// Pseudocode for line sweep
//
// Transform all coordinates (x, y) into (x + ιy, y) for infinitesimal ι to remove vertical lines and incidental same-time events
// Delete all length-0 segments
// PQ <- new priority queue
// SL <- new AVL tree
// add start and end events of all segments into PQ
// while e <- pop PQ
//     E <- pop all events happening at the same time as e from PQ
//     add e to E
//     SL_E <- remove all segments that are at the same position at e's time from SL
//     for all End(s) in E
//         remove s from SL_E
//     reverse order of SL_E
//     split segments in SL_E
//     for all Start(s) in E
//         insert s into SL_E
//     insert SL_E into SL
//     robustly check intersection between start of SL_E and previous in SL, and add Intersection event to PQ if there's one
//     if SL_E is not empty
//         robustly check intersection between end of SL_E and next in SL, and add Intersection event to PQ if there's one
// Now do it all again, but transform coordinates as (y + ιx, x) instead. (skip for exact coordinates)

#[derive(Clone, Debug)]
struct EventCmp<U, E>([U; 2], EventCase<E>);

impl<U: PartialEq, E> PartialEq for EventCmp<U, E> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<U: Eq, E> Eq for EventCmp<U, E> {}

impl<U: PartialOrd, E> PartialOrd for EventCmp<U, E> {
    /// Reversed for use in BinaryHeap
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        other.0.partial_cmp(&self.0)
    }
}

impl<U: Ord, E> Ord for EventCmp<U, E> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.0.cmp(&self.0)
    }
}

/// An `Option` with its comparison between Some and None reversed.
/// Needed because vertical lines have the highest slope after the (x + εy, y) transform.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum ReverseOption<T> {
    Some(T),
    None
}

/// I'm giving up and using a Vec here because there seems to be no Rust library
/// that supports a search tree with both
/// * the ability to use a comparator with mutable runtime information
/// * the ability to get the prev/next entry given a position
#[derive(Clone, Debug)]
struct SegmentSearchRef<E, T, F> {
    phantom: PhantomData<T>,
    segments: Vec<E>,
    mapping: F
}

impl<'a, E, T, F> SegmentSearchRef<E, T, F> where
    E: Eq + Clone,
    T: IntersectCoordinate<E>,
    F: Fn(&E) -> [VectorView2Dyn<'a, T>; 2] + Clone,
    for<'b> &'b T: RefNum<T>
{
    /// Tries to yoink a segment with the same position (not necessarily angle) as `s`.
    fn yoink_by_pos(&mut self, s: &E, curr_pos: <T as IntersectCoordinate<E>>::Time<'_>) -> Option<E> {
        let key = T::pos(s, self.mapping.clone(), curr_pos.clone());
        match self.segments.binary_search_by_key(&key, |seg| T::pos(seg, self.mapping.clone(), curr_pos.clone())) {
            Ok(pos) => Some(self.segments.remove(pos)),
            Err(_) => None,
        }
    }

    /// Inserts a segment by just angle.
    /// Presumably, the segments all tie by position
    /// 
    /// Returns `None` if succeeded.
    /// Returns the offending segment if there was already one there.
    fn insert_by_angle(&mut self, s: E) -> Option<E> {
        let key = T::angle(&s, self.mapping.clone());
        match self.segments.binary_search_by_key(&key, |seg| T::angle(seg, self.mapping.clone())) {
            Ok(pos) => Some(self.segments[pos].clone()),
            Err(pos) => { self.segments.insert(pos, s); None },
        }
    }

    /// Removes a segment
    fn remove(&mut self, s: &E) {
        // Just linear search for now, because this is called only on SegmentSearchRef with all segments at the same position
        // (though not necessarily angle)
        if let Some(pos) = self.segments.iter().position(|seg| seg == s) {
            self.segments.remove(pos);
        }
    }

    fn segments(&self) -> &[E] {
        &self.segments
    }

    fn reverse(&mut self) {
        self.segments.reverse();
    }

    /// In the result, all the segments are tied for position at `time`
    fn separate_all_with_same_pos(&mut self, s: &E, time: <T as IntersectCoordinate<E>>::Time<'_>) -> (Self, usize) {
        let key = T::pos(&s, self.mapping.clone(), time.clone());
        match self.segments.binary_search_by_key(&key, |seg| T::pos(seg, self.mapping.clone(), time.clone())) {
            Ok(pos) => {
                let mut lower_bound = pos;
                let mut upper_bound = pos + 1;
                while self.segments.get(lower_bound - 1).map(|s| T::pos(s, self.mapping.clone(), time.clone()) == key).unwrap_or(false) {
                    lower_bound -= 1;
                }
                while self.segments.get(upper_bound).map(|s| T::pos(s, self.mapping.clone(), time.clone()) == key).unwrap_or(false) {
                    upper_bound += 1;
                }
                let segments = self.segments.drain(lower_bound..upper_bound).collect::<Vec<_>>();
                (Self { phantom: PhantomData, segments, mapping: self.mapping.clone() }, lower_bound)
            }, Err(pos) => (Self { phantom: PhantomData, segments: vec![], mapping: self.mapping.clone() }, pos)
        }
    }

    /// Returns the end position of the extension (the start position is an argument; you already know that)
    fn extend(&mut self, other: Self, pos: usize) -> usize {
        let end = self.segments.split_off(pos);
        self.segments.extend(other.segments);
        let len = self.segments.len();
        self.segments.extend(end);
        len
    }

    /// Previous segment at position
    fn prev(&self, pos: usize) -> Option<&E> {
        if pos == 0 { return None }
        self.segments.get(pos - 1)
    }

    /// Next segment at position
    fn next(&self, pos: usize) -> Option<&E> {
        self.segments.get(pos)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EventCase<E> {
    Start(E),
    Intersect([E; 2]),
    End(E),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Endpoint {
    Start,
    End,
}

#[derive(Clone, Copy, Debug)]
pub enum EventF64<E> {
    Endpoint(Vector2<f64>, Endpoint, E),
    Intersect([[Vector2<f64>; 2]; 2], [E; 2]),
}

impl<E: PartialEq> PartialEq for EventF64<E> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Endpoint(t1, _, _), Self::Endpoint(t2, _, _)) => t1 == t2,
            (Self::Endpoint(t, _, _), Self::Intersect(segs, _)) => { todo!() },
            (Self::Intersect(segs, _), Self::Endpoint(t, _, _)) => other == self,
            (Self::Intersect(segs1, _), Self::Intersect(segs2, _)) => { todo!() },
        }
    }
}

impl<E: Eq> Eq for EventF64<E> {}

impl<E: PartialOrd> PartialOrd for EventF64<E> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self, other) {
            (Self::Endpoint(t1, _, _), Self::Endpoint(t2, _, _)) => {
                let t1 = [FloatOrd(t1.x), FloatOrd(t1.y)];
                let t2 = [FloatOrd(t2.x), FloatOrd(t2.y)];
                Some(t1.cmp(&t2))
            },
            (Self::Endpoint(t, _, _), Self::Intersect(segs, _)) => { todo!() },
            (Self::Intersect(segs, _), Self::Endpoint(t, _, _)) => other.partial_cmp(self).map(|o| o.reverse()),
            (Self::Intersect(segs1, _), Self::Intersect(segs2, _)) => { todo!() },
        }
    }
}

impl<E: Ord> Ord for EventF64<E> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[derive(Clone, Debug)]
pub enum EventExact<'a, E> {
    Endpoint(VectorView2Dyn<'a, BasedExpr>, Endpoint, E),
    Intersect(Vector2<BasedExpr>, [E; 2]),
}

impl<'a, E> EventExact<'a, E> {
    pub fn time(&'a self) -> VectorView2Dyn<'a, BasedExpr> {
        match self {
            Self::Endpoint(t, _, _) => t.clone(),
            Self::Intersect(t, _) => t.as_view()
        }
    }
}

impl<'a, E: PartialEq> PartialEq for EventExact<'a, E> {
    fn eq(&self, other: &Self) -> bool {
        self.time() == other.time()
    }
}

impl<'a, E: Eq> Eq for EventExact<'a, E> {}

impl<'a, E: PartialOrd> PartialOrd for EventExact<'a, E> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Safety: this is a VectorView2 ([2, 1]-matrix) with row stride 1,
        // so the elements are contiguous.
        Some(unsafe { self.time().data.as_slice_unchecked().cmp(other.time().data.as_slice_unchecked()) })
    }
}

impl<'a, E: Ord> Ord for EventExact<'a, E> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

pub trait Event: Sized {
    type This<'a>: Ord + Event<E = Self::E, N = Self::N>;
    type E: Clone + Ord;
    type N: IntersectCoordinate<Self::E>;

    fn try_start_end<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(segment: Self::E, mapping: F) -> Option<[Self::This<'a>; 2]>;
    fn try_intersect<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(segments: [Self::E; 2], mapping: F) -> Option<Self::This<'a>>;
    fn involved_segments(&self) -> &[Self::E];
    fn time<'a>(&'a self) -> <Self::N as IntersectCoordinate<Self::E>>::Time<'a>;
    fn segment_if_end(&self) -> Option<&Self::E>;
    fn segment_if_start(&self) -> Option<&Self::E>;
}

impl<E: Clone + Ord> Event for EventF64<E> {
    type This<'a> = Self;
    type E = E;
    type N = f64;

    fn try_start_end<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(segment: Self::E, mut mapping: F) -> Option<[Self::This<'a>; 2]> {
        let [mut start, mut end] = mapping(&segment);
        match [FloatOrd(start.x), FloatOrd(start.y)].cmp(&[FloatOrd(end.x), FloatOrd(end.y)]) {
            Ordering::Less => {},
            Ordering::Equal => return None, // length is 0
            Ordering::Greater => mem::swap(&mut start, &mut end),
        }
        Some([
            Self::Endpoint(start.into_owned(), Endpoint::Start, segment.clone()),
            Self::Endpoint(end.into_owned(), Endpoint::End, segment)
        ])
    }

    fn try_intersect<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(segments: [Self::E; 2], mut mapping: F) -> Option<Self::This<'a>> {
        let lines = segments.clone().map(|s| mapping(&s));
        todo!(); // TODO: Robust predicate Orient2D
        Some(Self::Intersect(lines.map(|line| line.map(|p| p.into_owned())), segments))
    }

    fn involved_segments(&self) -> &[Self::E] {
        match self {
            Self::Endpoint(_, _, e) => slice::from_ref(e),
            Self::Intersect(_, es) => es
        }
    }

    fn time<'a>(&'a self) -> <Self::N as IntersectCoordinate<Self::E>>::Time<'a> {
        match self {
            Self::Endpoint(t , _, _) => TimeF64::Point(*t),
            Self::Intersect(t, _) => TimeF64::Intersect(*t)
        }
    }

    fn segment_if_end(&self) -> Option<&Self::E> {
        if let Self::Endpoint(_, Endpoint::End, e) = self { Some(e) } else { None }
    }

    fn segment_if_start(&self) -> Option<&Self::E> {
        if let Self::Endpoint(_, Endpoint::Start, e) = self { Some(e) } else { None }
    }
}

impl<'b, E: Clone + Ord> Event for EventExact<'b, E> {
    type This<'a> = EventExact<'a, E>;
    type E = E;
    type N = BasedExpr;

    fn try_start_end<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(segment: Self::E, mut mapping: F) -> Option<[Self::This<'a>; 2]> {
        let [mut start, mut end] = mapping(&segment);
        // Safety: this is a VectorView2 ([2, 1]-matrix) with row stride 1,
        // so the elements are contiguous.
        match unsafe { start.data.as_slice_unchecked().cmp(end.data.as_slice_unchecked()) } {
            Ordering::Less => {},
            Ordering::Equal => return None, // length is 0
            Ordering::Greater => mem::swap(&mut start, &mut end),
        }
        Some([<Self::This<'a>>::Endpoint(start, Endpoint::Start, segment.clone()), <Self::This<'a>>::Endpoint(end, Endpoint::End, segment.clone())])
    }

    fn try_intersect<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(segments: [Self::E; 2], mut mapping: F) -> Option<Self::This<'a>> {
        let [line_a, line_b] = segments.clone().map(|s| mapping(&s));
        if let SegmentIntersection::Intersection(point) = geom::segment_intersect(line_a, line_b) {
            Some(<Self::This<'a>>::Intersect(point, segments))
        } else {
            None
        }
    }

    fn involved_segments(&self) -> &[Self::E] {
        match self {
            Self::Endpoint(_, _, e) => slice::from_ref(e),
            Self::Intersect(_, es) => es
        }
    }

    fn time<'a>(&'a self) -> <Self::N as IntersectCoordinate<Self::E>>::Time<'a> {
        match self {
            Self::Endpoint(t, _, _) => t.clone(),
            Self::Intersect(t, _) => t.as_view(),
        }
    }

    fn segment_if_end(&self) -> Option<&Self::E> {
        if let Self::Endpoint(_, Endpoint::End, e) = self { Some(e) } else { None }
    }

    fn segment_if_start(&self) -> Option<&Self::E> {
        if let Self::Endpoint(_, Endpoint::Start, e) = self { Some(e) } else { None }
    }
}

pub trait Time: Clone {
    type N: NumEx + RealField;
    /// Returns an Option even though the position should exist.
    /// 
    /// In the rare case that two lines intersect according to the robust geometric predicate,
    /// but the point can't be found due to a matrix that ends up being singular by say,
    /// catastrophic cancellation, we don't want to panic.
    /// 
    /// Note that the return value is always `Some(_)` for exact coordinates.
    fn actual_pos(self) -> Option<Vector2<Self::N>>;
}

impl Time for TimeF64 {
    type N = f64;
    fn actual_pos(self) -> Option<Vector2<Self::N>> {
        match self {
            Self::Point(p) => Some(p),
            Self::Intersect([[s1p0, s1p1], [s2p0, s2p1]]) => {
                if let LineIntersection::Intersection(p) = geom::line_intersect(
                    [s1p0.as_view(), s1p1.as_view()],
                    [s2p0.as_view(), s2p1.as_view()],
                ) { Some(p) } else { None }
            }
        }
    }
}

impl<'a> Time for VectorView2Dyn<'a, BasedExpr> {
    type N = BasedExpr;
    
    fn actual_pos(self) -> Option<Vector2<Self::N>> {
        Some(self.into_owned())
    }
}

#[derive(Clone, Copy, Debug)]
pub enum TimeF64 {
    Point(Vector2<f64>),
    Intersect([[Vector2<f64>; 2]; 2])
}

#[derive(Clone, Copy, Debug)]
pub struct PosF64([Vector2<f64>; 2], TimeF64);

impl PartialEq for PosF64 {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for PosF64 {}

impl PartialOrd for PosF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PosF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        todo!()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct AngleF64([Vector2<f64>; 2]);

impl PartialEq for AngleF64 {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for AngleF64 {}

impl PartialOrd for AngleF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AngleF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        todo!()
    }
}

pub trait IntersectCoordinate<E>: NumEx + RealField {
    type Event<'a>: Ord + Event<E = E, N = Self>;
    type Pos: Ord;
    type Time<'a>: Time<N = Self>;
    type Angle: Ord;

    /// The position of `segment` where the sweep line is at `curr_pos`
    fn pos<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mapping: F, time: Self::Time<'_>) -> Self::Pos;
    /// The angle of `segment`, pointing right of the sweep line
    fn angle<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mapping: F) -> Self::Angle;
}

impl<E: Ord + Clone> IntersectCoordinate<E> for f64 {
    type Event<'a> = EventF64<E>;
    type Pos = PosF64;
    type Time<'a> = TimeF64;
    type Angle = AngleF64;

    fn pos<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mut mapping: F, time: Self::Time<'_>) -> Self::Pos {
        let points = mapping(segment);
        PosF64(points.map(|p| p.into_owned()), time)
    }

    fn angle<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mut mapping: F) -> Self::Angle {
        let [mut p0, mut p1] = mapping(segment);
        if [FloatOrd(p0.x), FloatOrd(p0.y)] > [FloatOrd(p1.x), FloatOrd(p1.y)] {
            mem::swap(&mut p0, &mut p1);
        }
        AngleF64([p0.into_owned(), p1.into_owned()])
    }
}

impl<E: Ord + Clone> IntersectCoordinate<E> for BasedExpr {
    type Event<'a> = EventExact<'a, E>;
    type Pos = [BasedExpr; 2];
    type Time<'a> = VectorView2Dyn<'a, BasedExpr>;
    type Angle = ReverseOption<BasedExpr>;

    fn pos<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mut mapping: F, time: Self::Time<'_>) -> Self::Pos {
        let [p0, p1] = mapping(segment);
        let [[x, y]] = if p0.x == p1.x {
            time.into_owned()
        } else {
            let t = (&time.x - &p0.x) / (&p1.x - &p0.x);
            (&p1 - &p0) * t + p0
        }.data.0;
        [x.into(), y.into()]
    }

    fn angle<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mut mapping: F) -> Self::Angle {
        let [p0, p1] = mapping(segment);
        let denom = &p1.x - &p0.x;
        if denom.is_zero() {
            ReverseOption::None
        } else {
            ReverseOption::Some(((&p1.y - &p0.y) / denom).into())
        }
    }
}

impl<E> EventCase<E> {
    fn involved_segments(&self) -> &[E] {
        match self {
            Self::Start(s) => slice::from_ref(s),
            Self::Intersect(ss) => ss,
            Self::End(s) => slice::from_ref(s),
        }
    }
}

pub fn intersect_all_segments_ref<'a,
    E: Eq + Clone + 'a,
    T: IntersectCoordinate<E>,
    F: Fn(&E) -> [VectorView2Dyn<'a, T>; 2] + Clone
>(edges: impl IntoIterator<Item = E>, mapping: F) -> Vec<(E, Vector2<T>)> where
    for<'b> &'b T: RefNum<T>
{
    let mut events = BinaryHeap::new();
    let mut segments = SegmentSearchRef {
        phantom: PhantomData,
        segments: vec![],
        mapping: mapping.clone()
    };
    let mut splits = vec![];
    // Add all start/end events
    for edge in edges {
        T::Event::try_start_end(edge, mapping.clone()).map(|evs| events.extend(evs));
    }

    while let Some(event) = events.pop() {
        // get all events happening at the same time
        let mut curr_events = vec![event];
        loop {
            let next = if let Some(ev) = events.peek() { ev } else { break };
            if next == &curr_events[0] {
                curr_events.push(events.pop().unwrap())
            } else { break }
        }
        let segment = &curr_events[0].involved_segments()[0];
        let time = curr_events[0].time();
        // This is guaranteed to get all the segments with events that happen at `time`,
        // as well as yoink segments that need to be yoinked according to the pseudocode
        let (mut curr_segments, insert_pos) = segments.separate_all_with_same_pos(segment, time.clone());

        // Remove all segments with End events
        for event in &curr_events {
            if let Some(s) = event.segment_if_end() {
                curr_segments.remove(s);
            }
        }

        // Reverse & split the rest
        curr_segments.reverse();
        if let Some(time_pos) = time.actual_pos() {
            for s in curr_segments.segments() {
                splits.push((s.clone(), time_pos.clone()));
            }
        }

        // Insert all segments with Start events
        for event in curr_events {
            if let Some(s) = event.segment_if_start() {
                curr_segments.insert_by_angle(s.clone());
            }
        }

        // Recombine and check for new events
        let upper_pos = segments.extend(curr_segments, insert_pos);
        segments.prev(insert_pos).zip(segments.next(insert_pos))
            .and_then(|(s0, s1)| T::Event::try_intersect([s0.clone(), s1.clone()], mapping.clone()))
            .map(|ev| events.push(ev));
        if insert_pos != upper_pos {
            segments.prev(upper_pos).zip(segments.next(upper_pos))
                .and_then(|(s0, s1)| T::Event::try_intersect([s0.clone(), s1.clone()], mapping.clone()))
                .map(|ev| events.push(ev));
        }
    }
    
    splits
}

impl Frame {
    /// Adds vertices at all points of intersections between the edges of `self`.
    /// 
    /// This requires coordinates to be 2D, and will remove all face information. Beware. 
    /// 
    /// In particular:
    /// * Duplicate edges are merged. Even ones that showed up during processing.
    /// * If two vertices have the same coordinates, they get merged. Even ones that showed up during processing.
    /// * If two edges are collinear and intersect in their interiors, each edge is split
    ///     wherever it intersects a boundary point of the other edge. Then duplicate edges are merged.
    /// * If an edge intersects a boundary point of another edge, the first edge is split at that boundary point.
    /// * If two non-collinear edges intersect in their interiors, they're both split at the point of intersection.
    pub fn intersect_all_edges_generic<T: NumEx + Coordinate>(&mut self) {
    }
}