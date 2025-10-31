use std::{cmp::{Ordering, Reverse}, collections::BinaryHeap, marker::PhantomData, mem, slice, hash::Hash};

use exact_number::BasedExpr;
use indexmap::IndexMap;
use nalgebra::{vector, RawStorage, RealField, Vector2};
use num_traits::{RefNum, Zero};
use robust_geometry as robust;

use crate::{Edge, Frame, filter::Coordinate, geom::{self, FloatOrd, LineIntersection, NumEx, SegmentIntersection, VectorView2Dyn}};

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
//     robustly check intersection between start of SL_E and previous in SL, and add Intersection event to PQ if there's one in the future
//     if SL_E is not empty
//         robustly check intersection between end of SL_E and next in SL, and add Intersection event to PQ if there's one in the future
// Now do it all again, but transform coordinates as (y + ιx, x) instead. (skip for exact coordinates)

fn sort_segment(mut s: [Vector2<f64>; 2]) -> [Vector2<f64>; 2] {
    if [FloatOrd(s[0].x), FloatOrd(s[0].y)] > [FloatOrd(s[1].x), FloatOrd(s[1].y)] {
        s.swap(0, 1);
    }
    s
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
    /// Inserts a segment by just angle.
    /// Presumably, the segments all tie by position
    /// 
    /// Returns `None` if succeeded.
    /// Returns the offending segment if there was already one there.
    fn insert_by_angle(&mut self, s: E) {
        let key = T::angle(&s, self.mapping.clone());
        match self.segments.binary_search_by_key(&key, |seg| T::angle(seg, self.mapping.clone())) {
            Ok(pos) => { self.segments.insert(pos, s) },
            Err(pos) => { self.segments.insert(pos, s) },
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
                while lower_bound > 0 &&
                    self.segments.get(lower_bound - 1).map(|s| T::pos(s, self.mapping.clone(), time.clone()) == key).unwrap_or(false) {
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

impl<E: Clone + Ord> PartialEq for EventF64<E> {
    fn eq(&self, other: &Self) -> bool {
        self.time() == other.time()
    }
}

impl<E: Clone + Ord> Eq for EventF64<E> {}

impl<E: Clone + Ord> PartialOrd for EventF64<E> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<E: Clone + Ord> Ord for EventF64<E> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.time().cmp(&other.time())
    }
}

#[derive(Clone, Debug)]
pub enum EventExact<'a, E> {
    Endpoint(VectorView2Dyn<'a, BasedExpr>, Endpoint, E),
    Intersect(Vector2<BasedExpr>, [E; 2]),
}

impl<'a, E: Clone + Ord> PartialEq for EventExact<'a, E> {
    fn eq(&self, other: &Self) -> bool {
        self.time() == other.time()
    }
}

impl<'a, E: Clone + Ord> Eq for EventExact<'a, E> {}

impl<'a, E: Clone + Ord> PartialOrd for EventExact<'a, E> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a, E: Clone + Ord> Ord for EventExact<'a, E> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.time().cmp(other.time())
    }
}

pub trait Event: Sized {
    type This<'a>: Ord + Event<E = Self::E, N = Self::N>;
    type E: Clone + Ord;
    type N: IntersectCoordinate<Self::E>;

    fn try_start_end<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(segment: Self::E, mapping: F) -> Option<[Self::This<'a>; 2]>;
    /// Fails if there's no intersection or it's in the past/present
    fn try_intersect<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>
        (segments: [Self::E; 2], mapping: F, time: <Self::N as IntersectCoordinate<Self::E>>::Time<'_>) -> Option<Self::This<'a>>;
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

    fn try_intersect<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(
        segments: [Self::E; 2], mut mapping: F, time: <Self::N as IntersectCoordinate<Self::E>>::Time<'_>
    ) -> Option<Self::This<'a>> {
        let mut lines = segments.clone().map(|s| mapping(&s)).map(|line| sort_segment(line.map(|p| p.into_owned())));
        if robust::cross_2d(lines[1][0], lines[1][1], lines[0][0], lines[0][1]) < 0.0 {
            lines.swap(0, 1); // sort segments clockwise
        }
        if robust::orient_2d(lines[1][0], lines[1][1], lines[0][0]) >= 0.0 ||
            robust::orient_2d(lines[1][0], lines[1][1], lines[0][1]) <= 0.0 {
            return None; // does not intersect
        }
        if TimeF64::Intersect(lines) > time {
            Some(Self::Intersect(lines, segments))
        } else { None }
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

    fn try_intersect<'a, F: FnMut(&Self::E) -> [VectorView2Dyn<'a, Self::N>; 2]>(
        segments: [Self::E; 2], mut mapping: F, time: <Self::N as IntersectCoordinate<Self::E>>::Time<'_>
    ) -> Option<Self::This<'a>> {
        let [line_a, line_b] = segments.clone().map(|s| mapping(&s));
        if let SegmentIntersection::Intersection(point) = geom::segment_intersect(line_a, line_b) {
            if &point.data.0[0] > time { Some(<Self::This<'a>>::Intersect(point, segments)) } else { None }
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
            // Safety: VectorView2Dyn is a length-2 column vector with a row stride of 1,
            // so there are indeed 2 elements and they are indeed contiguous.
            Self::Endpoint(t, _, _) => t.as_slice().try_into().unwrap(),
            Self::Intersect(t, _) => &t.data.0[0],
        }
    }

    fn segment_if_end(&self) -> Option<&Self::E> {
        if let Self::Endpoint(_, Endpoint::End, e) = self { Some(e) } else { None }
    }

    fn segment_if_start(&self) -> Option<&Self::E> {
        if let Self::Endpoint(_, Endpoint::Start, e) = self { Some(e) } else { None }
    }
}

pub trait Time: Clone + Ord {
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

impl<'a> Time for &'a [BasedExpr; 2] {
    type N = BasedExpr;
    
    fn actual_pos(self) -> Option<Vector2<Self::N>> {
        Some(vector![self[0].clone(), self[1].clone()])
    }
}

/// For Intersect, stores left point before right point,
/// and stores segments in clockwise order
#[derive(Clone, Copy, Debug)]
pub enum TimeF64 {
    Point(Vector2<f64>),
    Intersect([[Vector2<f64>; 2]; 2])
}

impl PartialEq for TimeF64 {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for TimeF64 {}

impl PartialOrd for TimeF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TimeF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Point(p1), Self::Point(p2)) => [FloatOrd(p1.x), FloatOrd(p1.y)].cmp(&[FloatOrd(p2.x), FloatOrd(p2.y)]),
            (Self::Point(p1), Self::Intersect(ss)) =>
                FloatOrd(robust::simple_intersect_compare_i(*p1, ss[0][0], ss[0][1], ss[1][0], ss[1][1], 0)).cmp(&FloatOrd(0.0)).then_with(||
                FloatOrd(robust::simple_intersect_compare_i(*p1, ss[0][0], ss[0][1], ss[1][0], ss[1][1], 1)).cmp(&FloatOrd(0.0))),
            (Self::Intersect(_), Self::Point(_)) => other.cmp(self).reverse(),
            (Self::Intersect(x1), Self::Intersect(x2)) =>
                FloatOrd(robust::complex_intersect_compare_i(x1[0][0], x1[0][1], x1[1][0], x1[1][1], x2[0][0], x2[0][1], x2[1][0], x2[1][1], 0)).cmp(&FloatOrd(0.0)).then_with(||
                FloatOrd(robust::complex_intersect_compare_i(x1[0][0], x1[0][1], x1[1][0], x1[1][1], x2[0][0], x2[0][1], x2[1][0], x2[1][1], 1)).cmp(&FloatOrd(0.0)))
        }
    }
}

/// Stores left point before right point
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
        debug_assert_eq!(self.1, other.1); // never compare positions across different times
        //println!("Comparing positions");
        //crate::my_dbg!(inline &self);
        //crate::my_dbg!(inline other);
        // note: self.1 is always the current time,
        // so we don't have to worry about vertical comparison
        let (s1, s2) = (self.0, other.0);
        let orient = robust::cross_2d(s2[0], s2[1], s1[0], s1[1]);
        //crate::my_dbg!(inline orient);
        if orient == 0.0 {
            return FloatOrd(robust::orient_2d(s2[0], s2[1], s1[0])).cmp(&FloatOrd(0.0))
        }
        //if orient < 0.0 { mem::swap(&mut s1, &mut s2) } (turns out you need to do both this *and* flip the result, which cancel out)
        match &self.1 {
            TimeF64::Point(p) =>
                FloatOrd(robust::simple_intersect_compare_i(*p, s1[0], s1[1], s2[0], s2[1], 0)).cmp(&FloatOrd(0.0)).then_with(||
                FloatOrd(robust::simple_intersect_compare_i(*p, s1[0], s1[1], s2[0], s2[1], 1)).cmp(&FloatOrd(0.0))),
            TimeF64::Intersect(ss) =>
                FloatOrd(robust::complex_intersect_compare_i(ss[0][0], ss[0][1], ss[1][0], ss[1][1], s1[0], s1[1], s2[0], s2[1], 0)).cmp(&FloatOrd(0.0)).then_with(||
                FloatOrd(robust::complex_intersect_compare_i(ss[0][0], ss[0][1], ss[1][0], ss[1][1], s1[0], s1[1], s2[0], s2[1], 1)).cmp(&FloatOrd(0.0))),
        }
    }
}

/// Stores left point before right point
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
        let orient = robust::cross_2d(other.0[0], other.0[1], self.0[0], self.0[1]);
        FloatOrd(orient).cmp(&FloatOrd(0.0))
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
    
    /// Finds the splits between a bunch of segments in a way that's robust to perturbances
    fn intersect_all_segments_ref<'a,
        F: Fn(&E) -> [VectorView2Dyn<'a, Self>; 2] + Clone
    >(edges: impl IntoIterator<Item = E>, mapping: F, epsilon: &Self) -> Vec<(E, Vector2<Self>)> where
        E: Clone + Hash + 'a,
        for<'b> &'b Self: RefNum<Self>;
}

impl<E: Ord + Clone> IntersectCoordinate<E> for f64 {
    type Event<'a> = EventF64<E>;
    type Pos = PosF64;
    type Time<'a> = TimeF64;
    type Angle = AngleF64;

    fn pos<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mut mapping: F, time: Self::Time<'_>) -> Self::Pos {
        let points = sort_segment(mapping(segment).map(|p| p.into_owned()));
        PosF64(points, time)
    }

    fn angle<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mut mapping: F) -> Self::Angle {
        let seg = sort_segment(mapping(segment).map(|p| p.into_owned()));
        AngleF64(seg)
    }

    /// The tricky one.
    /// We add perpendicular serifs of length sqrt(2) * `epsilon` to each end of each segment,
    /// then extend each main segment (not the serifs) by sqrt(1/2) * `epsilon` on each side.
    /// Then, we remove all the serif splits.
    /// This could give erroneous splits; it is the job of the caller to merge splits close to each other.
    fn intersect_all_segments_ref<'a,
            F: Fn(&E) -> [VectorView2Dyn<'a, Self>; 2] + Clone
        >(edges: impl IntoIterator<Item = E>, mapping: F, epsilon: &Self) -> Vec<(E, Vector2<Self>)> where
            E: Clone + Hash + 'a,
            for<'b> &'b Self: RefNum<Self>
    {
        let edges = edges.into_iter().collect::<Vec<_>>();
        let serifs = edges.iter().map(|e| {
            let [p0, p1] = mapping(e);
            let diff = p1 - p0;
            let len = diff.norm();
            // also avoid dividing by really tiny numbers, like, idk, 5e-324.
            let unit = if len <= *epsilon { Vector2::zeros() } else { diff * *epsilon * (0.5f64).sqrt() / len };
            let perp = geom::perp_ccw(unit);
            (e.clone(), [p0 - perp, p0 - unit, p0 + perp, p1 + perp, p1 + unit, p1 - perp])
        }).collect::<IndexMap<_, _>>();

        let splits = intersect_all_segments_ref(
            edges.into_iter().flat_map(|e| (0..7).map(move |i| (e.clone(), i))),
            |(e, i)| {
                match i {
                    0 => mapping(e),
                    1..4 => [mapping(e)[0], serifs[e][i - 1].as_view()],
                    4..7 => [mapping(e)[1], serifs[e][i - 1].as_view()],
                    _ => unreachable!()
                }
            }
        );
        splits.into_iter()
            .filter(|((_, i), _) | *i == 0)
            .map(|((e, _), intersection)| (e, intersection))
            .collect()
    }
}

impl<E: Ord + Clone> IntersectCoordinate<E> for BasedExpr {
    type Event<'a> = EventExact<'a, E>;
    type Pos = [BasedExpr; 2];
    type Time<'a> = &'a [BasedExpr; 2];
    type Angle = ReverseOption<BasedExpr>;

    fn pos<'a, F: FnMut(&E) -> [VectorView2Dyn<'a, Self>; 2]>(segment: &E, mut mapping: F, time: Self::Time<'_>) -> Self::Pos {
        let [p0, p1] = mapping(segment);
        let [[x, y]] = if p0.x == p1.x {
            vector![time[0].clone(), time[1].clone()]
        } else {
            let t = (&time[0] - &p0.x) / (&p1.x - &p0.x);
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

    fn intersect_all_segments_ref<'a,
            F: Fn(&E) -> [VectorView2Dyn<'a, Self>; 2] + Clone
        >(edges: impl IntoIterator<Item = E>, mapping: F, _epsilon: &Self) -> Vec<(E, Vector2<Self>)> where
            E: Clone + Hash + 'a,
            for<'b> &'b Self: RefNum<Self>
    {
        intersect_all_segments_ref(edges, mapping)
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
        T::Event::try_start_end(edge, mapping.clone()).map(|evs| events.extend(evs.map(Reverse)));
    }

    while let Some(Reverse(event)) = events.pop() {
        //println!("");
        // get all events happening at the same time
        let mut curr_events = vec![event];
        loop {
            let next = if let Some(Reverse(ev)) = events.peek() { ev } else { break };
            if next == &curr_events[0] {
                curr_events.push(events.pop().unwrap().0)
            } else { break }
        }
        let segment = &curr_events[0].involved_segments()[0];
        let time = curr_events[0].time();
        // This is guaranteed to get all the segments with events that happen at `time`,
        // as well as yoink segments that need to be yoinked according to the pseudocode
        let (mut curr_segments, insert_pos) = segments.separate_all_with_same_pos(segment, time.clone());
        //crate::my_dbg!(inline &time);

        // Remove all segments with End events
        //crate::my_dbg!(inline time.clone().actual_pos());
        for event in &curr_events {
            //crate::my_dbg!(inline event);
            if let Some(s) = event.segment_if_end() {
                curr_segments.remove(s);
            }
        }

        //println!("Current segments");
        //for s in curr_segments.segments() {
        //    crate::my_dbg!(inline s);
        //}

        // Reverse & split the rest
        curr_segments.reverse();
        if let Some(time_pos) = time.clone().actual_pos() {
            for s in curr_segments.segments() {
                splits.push((s.clone(), time_pos.clone()));
            }
        }

        // Insert all segments with Start events
        for event in &curr_events {
            if let Some(s) = event.segment_if_start() {
                curr_segments.insert_by_angle(s.clone());
            }
        }

        // Recombine and check for new events
        let upper_pos = segments.extend(curr_segments, insert_pos);
        segments.prev(insert_pos).zip(segments.next(insert_pos))
            .and_then(|(s0, s1)| T::Event::try_intersect([s0.clone(), s1.clone()], mapping.clone(), time.clone()))
            .map(|ev| events.push(Reverse(ev)));
        if insert_pos != upper_pos {
            segments.prev(upper_pos).zip(segments.next(upper_pos))
                .and_then(|(s0, s1)| T::Event::try_intersect([s0.clone(), s1.clone()], mapping.clone(), time))
                .map(|ev| events.push(Reverse(ev)));
        }

        //println!("All segments");
        //for s in segments.segments() {
        //    crate::my_dbg!(inline s);
        //}
    }
    
    splits
}

impl Frame {
    /// Adds vertices at all points of intersections between the edges of `self`.
    /// 
    /// This requires coordinates to be 2D and edge data to exist and will remove all face information. Beware. 
    /// 
    /// In particular:
    /// * Duplicate edges are merged. Even ones that showed up during processing.
    /// * If two vertices have the same coordinates, they get merged. Even ones that showed up during processing.
    /// * If two edges are collinear and intersect in their interiors, each edge is split
    ///     wherever it intersects a boundary point of the other edge. Then duplicate edges are merged.
    /// * If an edge intersects a boundary point of another edge, the first edge is split at that boundary point.
    /// * If two non-collinear edges intersect in their interiors, they're both split at the point of intersection.
    pub fn intersect_all_edges_generic<T: NumEx + Coordinate + IntersectCoordinate<Edge>>(&mut self, epsilon: &T) where
        for<'a> &'a T: RefNum<T>
    {
        let vertices_coords = T::vertices_coords(self).as_ref().unwrap(); // required by spec
        assert_eq!(vertices_coords.nrows(), 2, "intersect_all_edges requires 2D coordinates");
        let edges_vertices = self.edges_vertices.as_ref().unwrap(); // required by spec
        let splits = T::intersect_all_segments_ref(
            (0..edges_vertices.len()).map(Edge),
            |&e| edges_vertices[e].map(|v| vertices_coords.fixed_view::<2, 1>(0, v.0)),
            epsilon
        );
    }
}

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use approx::relative_ne;
    use exact_number::{BasedExpr};
    use nalgebra::{vector, Vector2};

    use crate::{filter::intersect::intersect_all_segments_ref, geom::FloatOrd};

    macro_rules! exact_vec_2 {
        (($($a:tt)*), ($($b:tt)*)$(,)?) => {
            nalgebra::vector![exact_number::based_expr!($($a)*), exact_number::based_expr!($($b)*)]
        };
    }

    fn canonicalize_f64<E: Ord + Clone>(mut vec: Vec<(E, Vector2<f64>)>) -> Vec<(E, Vector2<f64>)> {
        fn key<'a, E>(split: &'a (E, Vector2<f64>)) -> (&'a E, [FloatOrd<f64>; 2]) {
            (&split.0, [FloatOrd(split.1.x), FloatOrd(split.1.y)])
        }
        vec.sort_by(|a, b| key(a).cmp(&key(b)));
        vec
    }

    fn canonicalize_exact<E: Ord + Clone>(mut vec: Vec<(E, Vector2<BasedExpr>)>) -> Vec<(E, Vector2<BasedExpr>)> {
        fn key<'a, E>(split: &'a (E, Vector2<BasedExpr>)) -> (&'a E, &'a [BasedExpr]) {
            (&split.0, split.1.data.as_slice())
        }
        vec.sort_by(|a, b| key(a).cmp(&key(b)));
        vec
    }

    fn round_vectors(vectors: Vec<Vector2<BasedExpr>>) -> Vec<Vector2<f64>> {
        vectors.into_iter()
            .map(|v| v.map(|c| c.round_to_nearest_f64()))
            .collect::<Vec<_>>()
    }

    fn round_splits<E>(splits: Vec<(E, Vector2<BasedExpr>)>) -> Vec<(E, Vector2<f64>)> {
        splits.into_iter()
            .map(|(e, v)| (e, v.map(|c| c.round_to_nearest_f64())))
            .collect::<Vec<_>>()
    }

    fn assert_splits<E: Debug + Eq>(splits: &Vec<(E, Vector2<f64>)>, expected: &Vec<(E, Vector2<f64>)>) {
        if splits.len() != expected.len() || splits.iter().zip(expected.iter())
            .any(|((s1e, s1v), (s2e, s2v))|
                s1e != s2e || relative_ne!(s1v, s2v))
        {
            panic!("splits assertion failed: {splits:?} is not close enough to {expected:?}")
        }
    }

    #[test]
    fn test_intersect_all_segments_one_dodge() {
        let vectors = vec![
            exact_vec_2![(0), (0)],
            exact_vec_2![(1), (1)],
            exact_vec_2![(0), (2)],
            exact_vec_2![(2), (1)],
        ];
        let segments = vec![
            [0, 1],
            [2, 3],
        ];
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![]);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![]);
    }

    #[test]
    fn test_intersect_all_segments_touches() {
        let vectors = vec![
            exact_vec_2![(0), (0)],
            exact_vec_2![(1), (1)],
            exact_vec_2![(0), (2)],
            exact_vec_2![(2), (0)],
            exact_vec_2![(3), (2)],
            exact_vec_2![(5), (0)],
            exact_vec_2![(5), (2)],
            exact_vec_2![(4), (1)],
        ];
        let segments = vec![
            [0, 1],
            [2, 3],
            [4, 5],
            [6, 7],
        ];
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (1, exact_vec_2![(1), (1)]),
            (2, exact_vec_2![(4), (1)]),
        ]);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (1, vector![1.0, 1.0]),
            (2, vector![4.0, 1.0]),
        ]);
    }

    #[test]
    fn test_intersect_all_segments_true_intersection() {
        let vectors = vec![
            exact_vec_2![(0), (0)],
            exact_vec_2![(1), (1)],
            exact_vec_2![(0), (1)],
            exact_vec_2![(2), (0)],
        ];
        let segments = vec![
            [0, 1],
            [2, 3],
        ];
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (0, exact_vec_2![(2/3), (2/3)]),
            (1, exact_vec_2![(2/3), (2/3)]),
        ]);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_splits(&splits, &vec![
            (0, vector![2.0 / 3.0, 2.0 / 3.0]),
            (1, vector![2.0 / 3.0, 2.0 / 3.0]),
        ]);
    }

    #[test]
    fn test_intersect_all_segments_multi_intersection() {
        let vectors = vec![
            exact_vec_2![(0), (0)],
            exact_vec_2![(1), (1)],
            exact_vec_2![(0), (1)],
            exact_vec_2![(1), (0)],
            exact_vec_2![(1/2), (1)],
            exact_vec_2![(1/2), (0)],
        ];
        let segments = vec![
            [0, 1],
            [2, 3],
            [4, 5],
        ];
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (0, exact_vec_2![(1/2), (1/2)]),
            (1, exact_vec_2![(1/2), (1/2)]),
            (2, exact_vec_2![(1/2), (1/2)]),
        ]);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (0, vector![0.5, 0.5]),
            (1, vector![0.5, 0.5]),
            (2, vector![0.5, 0.5]),
        ]);
    }

    #[test]
    fn test_intersect_all_segments_segment_stops_in_cove() {
        let vectors = vec![
            exact_vec_2![(0), (0)],
            exact_vec_2![(1), (1)],
            exact_vec_2![(0), (1)],
            exact_vec_2![(1), (0)],
            exact_vec_2![(1/3), (1/2)],
            exact_vec_2![(-1), (1/2)],
        ];
        let segments = vec![
            [0, 1],
            [2, 3],
            [4, 5], // this segment stops in the cove created by segments 0 and 1
        ];
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (0, exact_vec_2![(1/2), (1/2)]),
            (1, exact_vec_2![(1/2), (1/2)]),
        ]);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (0, vector![0.5, 0.5]),
            (1, vector![0.5, 0.5]),
        ]);
    }

    #[test]
    fn test_intersect_all_segments_segment_touches_cove() {
        let vectors = vec![
            exact_vec_2![(0), (0)],
            exact_vec_2![(1), (1)],
            exact_vec_2![(0), (1)],
            exact_vec_2![(1), (0)],
            exact_vec_2![(1/2), (1/2)],
            exact_vec_2![(-1), (1/2)],
            exact_vec_2![(3), (0)],
            exact_vec_2![(4), (1)],
            exact_vec_2![(3), (1)],
            exact_vec_2![(4), (0)],
            exact_vec_2![(7/2), (1/2)],
            exact_vec_2![(10), (1/2)],
        ];
        let segments = vec![
            [0, 1],
            [2, 3],
            [4, 5], // this segment touches the cove created by segments 0 and 1
            [6, 7],
            [8, 9],
            [10, 11],
        ];
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (0, exact_vec_2![(1/2), (1/2)]),
            (1, exact_vec_2![(1/2), (1/2)]),
            (3, exact_vec_2![(7/2), (1/2)]),
            (4, exact_vec_2![(7/2), (1/2)]),
        ]);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (0, vector![0.5, 0.5]),
            (1, vector![0.5, 0.5]),
            (3, vector![3.5, 0.5]),
            (4, vector![3.5, 0.5]),
        ]);
    }

    #[test]
    fn test_intersect_all_segments_coplanar() {
        let vectors = vec![
            exact_vec_2![(0), (0)],
            exact_vec_2![(0), (1)],
            exact_vec_2![(0), (0)],
            exact_vec_2![(0), (3)],
            exact_vec_2![(0), (2)],
            exact_vec_2![(0), (4)],
            exact_vec_2![(0), (4)],
            exact_vec_2![(0), (3)],
        ];
        let segments = vec![
            [0, 1],
            [2, 3],
            [4, 5],
            [6, 7],
        ];
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (1, exact_vec_2![(0), (1)]),
            (1, exact_vec_2![(0), (2)]),
            (2, exact_vec_2![(0), (3)]),
        ]);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits, vec![
            (1, vector![0.0, 1.0]),
            (1, vector![0.0, 2.0]),
            (2, vector![0.0, 3.0]),
        ]);
    }

    #[test]
    fn test_intersect_all_segments_big_intersection_predicates() {
        // Tests the predicates by making sure this returns the correct number of splits
        let vectors = vec![
            vector![0.0, 0.0],
            vector![1.0, 1.0],
            vector![0.0, 0.5],
            vector![1.0, 0.0],
            vector![0.5, 0.0],
            vector![0.0, 1.0],
            vector![0.25, 0.0],
            vector![0.5, 1.0],
            vector![0.0, 0.25],
            vector![1.0, 0.5],
            vector![-2.0, 0.0],
            vector![5.0, 1.0],
            vector![7.0, 6.0],
            vector![-13.0, -11.0],
            vector![123456789.0, 987654321.0],
            vector![-246913577.0, -1975308641.0],
        ];
        // All the segments intersect at (1/3, 1/3),
        // which can't be represented as an f64.
        // The algorithm should still tell that the intersection points are the same.
        // and make only 1 split per segment.
        let segments = vec![
            [0, 1],
            [2, 3],
            [4, 5],
            [6, 7],
            [8, 9],
            [10, 11],
            [12, 13],
            [14, 15]
        ];
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        assert_eq!(splits.len(), segments.len());
    }

    #[test]
    fn test_intersect_all_segments_small_grid() {
        let exact = |i: usize| BasedExpr::with_rational_basis(i.into());
        let size = 2;
        let vectors = (0..=size)
            .flat_map(|i| [vector![exact(0), exact(i)], vector![exact(size), exact(i)]])
            .chain((0..=size).flat_map( |i| [vector![exact(i), exact(0)], vector![exact(i), exact(size)]]))
            .collect::<Vec<_>>();
        let segments = (0..=size).map(|i| [2 * i, 2 * i + 1])
            .chain((0..=size).map(|i| [2 * (size + 1 + i), 2 * (size + 1 + i) + 1]))
            .collect::<Vec<_>>();
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        let expected = canonicalize_exact((0..=size).flat_map(|y| (0..=size).flat_map(move |x| {
            let mut result = vec![];
            if x > 0 && x < size {
                result.push((y, vector![exact(x), exact(y)]));
            }
            if y > 0 && y < size {
                result.push((size + 1 + x, vector![exact(x), exact(y)]))
            }
            result
        })).collect::<Vec<_>>());
        assert_eq!(splits, expected);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        let expected = round_splits(expected);
        assert_eq!(splits, expected);
    }

    #[test]
    fn test_intersect_all_segments_grid() {
        let exact = |i: usize| BasedExpr::with_rational_basis(i.into());
        let size = 8;
        let vectors = (0..=size)
            .flat_map(|i| [vector![exact(0), exact(i)], vector![exact(size), exact(i)]])
            .chain((0..=size).flat_map( |i| [vector![exact(i), exact(0)], vector![exact(i), exact(size)]]))
            .collect::<Vec<_>>();
        let segments = (0..=size).map(|i| [2 * i, 2 * i + 1])
            .chain((0..=size).map(|i| [2 * (size + 1 + i), 2 * (size + 1 + i) + 1]))
            .collect::<Vec<_>>();
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        let expected = canonicalize_exact((0..=size).flat_map(|y| (0..=size).flat_map(move |x| {
            let mut result = vec![];
            if x > 0 && x < size {
                result.push((y, vector![exact(x), exact(y)]));
            }
            if y > 0 && y < size {
                result.push((size + 1 + x, vector![exact(x), exact(y)]))
            }
            result
        })).collect::<Vec<_>>());
        assert_eq!(splits, expected);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        let expected = round_splits(expected);
        assert_eq!(splits, expected);
    }

    #[test]
    fn test_intersect_all_segments_hexagram() {
        // Now we're getting a little complicated
        let vectors = vec![
            exact_vec_2![(1 + 0 sqrt 3), (0 + 0 sqrt 3)],
            exact_vec_2![(1/2 + 0 sqrt 3), (0 + 1/2 sqrt 3)],
            exact_vec_2![(-1/2 + 0 sqrt 3), (0 + 1/2 sqrt 3)],
            exact_vec_2![(-1 + 0 sqrt 3), (0 + 0 sqrt 3)],
            exact_vec_2![(-1/2 + 0 sqrt 3), (0 - 1/2 sqrt 3)],
            exact_vec_2![(1/2 + 0 sqrt 3), (0 - 1/2 sqrt 3)],
        ];
        let segments = vec![
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
            [0, 2], [1, 3], [2, 4], [3, 5], [4, 0], [5, 1],
            [0, 3], [1, 4], [2, 5],
        ];
        let splits = canonicalize_exact(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        let expected = vec![
            (6, exact_vec_2![(0 + 0 sqrt 3), (0 + 1/3 sqrt 3)]),
            (6, exact_vec_2![(1/4 + 0 sqrt 3), (0 + 1/4 sqrt 3)]),
            (6, exact_vec_2![(1/2 + 0 sqrt 3), (0 + 1/6 sqrt 3)]),
            (7, exact_vec_2![(-1/2 + 0 sqrt 3), (0 + 1/6 sqrt 3)]),
            (7, exact_vec_2![(-1/4 + 0 sqrt 3), (0 + 1/4 sqrt 3)]),
            (7, exact_vec_2![(0 + 0 sqrt 3), (0 + 1/3 sqrt 3)]),
            (8, exact_vec_2![(-1/2 + 0 sqrt 3), (0 - 1/6 sqrt 3)]),
            (8, exact_vec_2![(-1/2 + 0 sqrt 3), (0 + 0 sqrt 3)]),
            (8, exact_vec_2![(-1/2 + 0 sqrt 3), (0 + 1/6 sqrt 3)]),
            (9, exact_vec_2![(-1/2 + 0 sqrt 3), (0 - 1/6 sqrt 3)]),
            (9, exact_vec_2![(-1/4 + 0 sqrt 3), (0 - 1/4 sqrt 3)]),
            (9, exact_vec_2![(0 + 0 sqrt 3), (0 - 1/3 sqrt 3)]),
            (10, exact_vec_2![(0 + 0 sqrt 3), (0 - 1/3 sqrt 3)]),
            (10, exact_vec_2![(1/4 + 0 sqrt 3), (0 - 1/4 sqrt 3)]),
            (10, exact_vec_2![(1/2 + 0 sqrt 3), (0 - 1/6 sqrt 3)]),
            (11, exact_vec_2![(1/2 + 0 sqrt 3), (0 - 1/6 sqrt 3)]),
            (11, exact_vec_2![(1/2 + 0 sqrt 3), (0 + 0 sqrt 3)]),
            (11, exact_vec_2![(1/2 + 0 sqrt 3), (0 + 1/6 sqrt 3)]),
            (12, exact_vec_2![(-1/2 + 0 sqrt 3), (0 + 0 sqrt 3)]),
            (12, exact_vec_2![(0 + 0 sqrt 3), (0 + 0 sqrt 3)]),
            (12, exact_vec_2![(1/2 + 0 sqrt 3), (0 + 0 sqrt 3)]),
            (13, exact_vec_2![(-1/4 + 0 sqrt 3), (0 - 1/4 sqrt 3)]),
            (13, exact_vec_2![(0 + 0 sqrt 3), (0 + 0 sqrt 3)]),
            (13, exact_vec_2![(1/4 + 0 sqrt 3), (0 + 1/4 sqrt 3)]),
            (14, exact_vec_2![(-1/4 + 0 sqrt 3), (0 + 1/4 sqrt 3)]),
            (14, exact_vec_2![(0 + 0 sqrt 3), (0 + 0 sqrt 3)]),
            (14, exact_vec_2![(1/4 + 0 sqrt 3), (0 - 1/4 sqrt 3)]),
        ];
        assert_eq!(splits, expected);

        let vectors = round_vectors(vectors);
        let splits = canonicalize_f64(intersect_all_segments_ref(0..segments.len(),
            |s| segments[*s].map(|v| vectors[v].as_view())));
        let expected = round_splits(expected);
        // The center intersection is the only one where >2 segments intersect,
        // and the way the coordinates are set up, (and the rounding mode of BasedExpr::round_to_nearest_f64),
        // those segments should intersect exactly at the same point,
        // so there are only 3 splits there, not 6.
        assert_splits(&splits, &expected);
    }
}