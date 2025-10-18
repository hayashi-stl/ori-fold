#![cfg(test)]

use typed_index_collections::{TiSlice, TiVec};

use crate::{fold::{AtFaceCorner, Edge, EdgesFaceCornersEx, Face, FaceCorner, Frame, HalfEdge, Vertex}, topology};

pub(crate) fn assert_ec_fh_consistent(ec: &TiSlice<Edge, [Vec<FaceCorner>; 2]>, fh: &TiSlice<Face, Vec<HalfEdge>>) {
    let rip = || panic!("{ec:?} is not consistent with {fh:?}");

    let mut fh_unaccounted = fh.to_vec();
    for (h, corners) in ec.half_iter_enumerated() {
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

impl Frame {
    /// Adds topology to this frame
    /// from data that guarantees a unique representation.
    /// `edges_face_corners` is sorted.
    /// without specifying it explicitly)
    pub(crate) fn with_topology_vh_fh(
        mut self,
        vertices_half_edges: Option<TiVec<Vertex, Vec<HalfEdge>>>,
        faces_half_edges: Option<TiVec<Face, Vec<HalfEdge>>>)
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
        self.vertices_half_edges = vertices_half_edges;
        self.edges_vertices = edges_vertices;
        self.edges_face_corners = edges_face_corners;
        self.faces_half_edges = faces_half_edges;
        self
    }
}