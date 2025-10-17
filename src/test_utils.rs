#![cfg(test)]

use typed_index_collections::TiVec;

use crate::{fold::{Face, Frame, HalfEdge, Vertex}, topology};

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