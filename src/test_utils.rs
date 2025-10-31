#![cfg(test)]

use indexmap::{indexset};
use typed_index_collections::{ti_vec};

use crate::{fold::{AtFaceCorner, Edge, EdgesFaceCornersEx, EdgesFaceCornersSlice, EdgesVerticesEx, EdgesVerticesSlice, FacesHalfEdges, FacesHalfEdgesSlice, Frame, FrameAttribute, HalfEdge, Vertex, VerticesHalfEdges, VerticesHalfEdgesSlice}, topology};

pub(crate) fn assert_ec_fh_consistent(ec: &EdgesFaceCornersSlice, fh: &FacesHalfEdgesSlice) {
    let rip = || panic!("{ec:?} is not consistent with {fh:?}");

    let mut fh_unaccounted = fh.to_vec();
    for (_, corners) in ec.half_iter_enumerated() {
        for &c in corners.iter() {
            if let Some(h) = fh_unaccounted.try_at_mut(c) {
                *h = HalfEdge(usize::MAX)
            } else {
                rip();
            }
        }
    }

    if fh_unaccounted.into_iter().any(|hs| hs.into_iter().any(|h| h != HalfEdge(usize::MAX))) {
        rip();
    }
}

pub(crate) fn assert_vh_ev_consistent(vh: &VerticesHalfEdgesSlice, ev: &EdgesVerticesSlice) {
    let rip = || panic!("{vh:?} is not consistent with {ev:?}");

    let mut ev_unaccounted = ev.to_vec();
    for half_edges in vh {
        for &h in half_edges.iter() {
            if let Some(v) = ev_unaccounted.try_at_mut(h) {
                *v[0] = Vertex(usize::MAX)
            } else {
                rip();
            }
        }
    }

    if ev_unaccounted.into_iter().any(|vs| vs.into_iter().any(|v| v != Vertex(usize::MAX))) {
        rip();
    }
}

impl Frame {
    pub(crate) fn assert_manifold_repr(&self) {
        let vh = self.vertices_half_edges.as_ref().unwrap();
        let ec = self.edges_face_corners.as_ref().unwrap();
        //let fh = self.faces_half_edges.as_ref().unwrap();
        for edge in (0..ec.len()).map(Edge) {
            if !self.is_edge_manifold(edge) {
                panic!("{self:?} is not in manifold representation: {edge} is not manifold");
            }
        }

        for (v, half_edges) in vh.iter_enumerated() {
            if let Some(fans) = self.vertex_fans_if_manifold(v) {
                let mut just_loop = false;
                let mut fan_map = ti_vec![indexset![]; ec.len() * 2];
                for (fan, is_loop) in fans {
                    just_loop = just_loop || is_loop;
                    for (&h0, &h1) in fan.iter().zip(fan.iter().cycle().skip(1)).take(fan.len() - if is_loop { 0 } else { 1 }) {
                        fan_map[h0].insert(h1);
                        fan_map[h1].insert(h0);
                    }
                }

                // Check that fans line up
                for (&h0, &h1) in half_edges.iter().zip(half_edges.iter().cycle().skip(1)).take(half_edges.len() - if just_loop { 0 } else { 1 }) {
                    fan_map[h0].swap_remove(&h1);
                    fan_map[h1].swap_remove(&h0);
                }
                if fan_map.into_iter().any(|set| !set.is_empty()) {
                    panic!("{self:?} is not in manifold representation: {v}'s half-edges are not in manifold order");
                }
            } else {
                panic!("{self:?} is not in manifold representation: {v} is not even manifold!");
            }
        }
    }

    pub(crate) fn assert_oriented(&self) {
        self.assert_manifold_repr();
        let vh = self.vertices_half_edges.as_ref().unwrap();
        let ec = self.edges_face_corners.as_ref().unwrap();
        let fh = self.faces_half_edges.as_ref().unwrap();
        for edge in (0..ec.len()).map(Edge) {
            if !self.is_edge_oriented(edge) {
                panic!("{self:?} is not oriented: {edge} is not oriented");
            }
        }

        let is_next = |h0: HalfEdge, h1: HalfEdge|
            ec.at(h0).len() > 0 && fh.at(ec.at(h0)[0].prev(fh)).flipped() == h1;

        for (v, half_edges) in vh.iter_enumerated() {
            // We don't need to check the first and last in case of a loop here.
            // If everything else is in order, then it just works.
            for (&h0, &h1) in half_edges.iter().zip(half_edges.iter().skip(1)) {
                if is_next(h1, h0) && !is_next(h0, h1) {
                    panic!("{self:?} is not oriented: {v} is not oriented");
                }
            }
        }
    }

    pub(crate) fn assert_topology_consistent(&self) {
        if let Some(edges_vertices) = self.edges_vertices.as_ref() {
            assert_eq!(self.num_vertices, self.vertices_half_edges.as_ref().unwrap().len());
            assert_vh_ev_consistent(self.vertices_half_edges.as_ref().unwrap(), &edges_vertices);
        }
        if let Some(faces_half_edges) = self.faces_half_edges.as_ref() {
            assert_ec_fh_consistent(self.edges_face_corners.as_ref().unwrap(), &faces_half_edges);
        }

        if self.faces_half_edges.is_some() {
            if self.frame_attributes.contains(&FrameAttribute::Orientable) {
                self.assert_oriented();
            } else if self.frame_attributes.contains(&FrameAttribute::Manifold) {
                self.assert_manifold_repr();
            }
        }
    }

    /// Adds topology to this frame
    /// from data that guarantees a unique representation.
    /// `edges_face_corners` is sorted.
    pub(crate) fn with_topology_vh_fh(
        mut self,
        vertices_half_edges: Option<VerticesHalfEdges>,
        faces_half_edges: Option<FacesHalfEdges>)
    -> Self {
        let edges_vertices = vertices_half_edges.as_ref().map(|vh| topology::try_vh_to_ev(&vh).unwrap());
        let mut edges_face_corners =
            faces_half_edges.as_ref().map(|fh|
                topology::try_fh_to_ec(&fh, edges_vertices.as_ref().unwrap().len()).unwrap());
        edges_face_corners.as_mut().map(|ec|
            ec.iter_mut().for_each(|cs| {
                cs[0].sort();
                cs[1].sort();
            }));
        if let Some(vh) = vertices_half_edges.as_ref() {
            self.num_vertices = vh.len();
        }
        self.vertices_half_edges = vertices_half_edges;
        self.edges_vertices = edges_vertices;
        self.edges_face_corners = edges_face_corners;
        self.faces_half_edges = faces_half_edges;
        self
    }
}

//fn setup_assert_oriented_prev_fail() -> Frame {
//    use FrameAttribute::*;
//    use HalfEdge as H;
//    use Vertex as V;
//    use Edge as E;
//
//    let mut frame = Frame {
//        faces_custom: indexmap! {
//            "test:test".to_owned() => ti_vec![json!(1)]
//        },
//        ..Default::default()
//    }.with_topology_vh_fh(Some(ti_vec![
//        vec![H(0), H(6)],
//        vec![H(7), H(2), H(8)],
//        vec![H(9), H(4)],
//        vec![H(10), H(1)],
//        vec![H(12), H(3), H(11)],
//        vec![H(5), H(13)],
//    ]), Some(ti_vec![
//        vec![H(2), H(12), H(5), H(9)],
//    ])).try_into_orientable().unwrap().0;
//
//    let index = frame.add_face(FaceData {
//        half_edges: vec![H(0), H(10), H(3), H(7)],
//        custom: indexmap! { "test:test".to_owned() => json!(3) }
//    });
//    let two_panel = frame;
//
//    let mut frame = two_panel.clone();
//    // Still orientable, but needs a flip
//    let [h14, h15] = frame.add_edge_like(E(0), [V(2), V(0)]).split();
//    let [h16, h17] = frame.add_edge_like(E(0), [V(5), V(3)]).split();
//    let index = frame.add_face(FaceData {
//        half_edges: vec![H(5), h14, H(0), h17],
//        custom: indexmap! { "test:test".to_owned() => json!(3) }
//    });
//    assert_eq!(frame.frame_attributes, indexset! {Manifold});
//    frame.vertices_half_edges.as_mut().unwrap()[V(0)].reverse();
//    frame
//}
//
//#[test]
//fn test_assert_oriented_prev_fail_early() {
//    setup_assert_oriented_prev_fail();
//}
//
//#[test]
//#[should_panic(expected = "is not oriented")]
//fn test_assert_oriented_prev_fail() {
//    let frame = setup_assert_oriented_prev_fail();
//    frame.assert_topology_consistent();
//}