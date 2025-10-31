use nalgebra::DVector;

use crate::{AtFaceCorner, Edge, EdgeData, EdgesFaceCornersEx, FaceCorner, Frame, FrameAttribute, HalfEdge, Vertex, filter::Coordinate, geom::NumEx};

impl Frame {
    /// Splits an edge at point `at`.
    /// If there are no vertex coordinates of type `T`, `at` gets ignored.
    /// The new vertex has the same data as vertex `edges_vertices[edge][0]`,
    /// except for coordinates. The new edge has the same data as `edge`
    /// except for the vertices it's connected to.
    /// 
    /// Returns the new vertex `v` created, and the new edge created from the split,
    /// which goes from `v` to `edges_vertices[edge][1]`.
    /// 
    /// This is a purely topological operation and might not preserve
    /// geometric properties such as planarity if you're not careful.
    /// It does preserve manifold-ness and oriented-ness.
    /// 
    /// Requires edge data to exist.
    pub fn split_edge<T: NumEx + Coordinate>(&mut self, edge: Edge, at: DVector<T>) -> (Vertex, Edge) {
        // Add new vertex
        let new_v = self.add_vertex_like(self.edges_vertices.as_ref().unwrap()[edge][0], at);

        // Find where the edge is stored in the end vertex, so it can be replaced later
        let ev = self.edges_vertices.as_ref().unwrap();
        let vh = self.vertices_half_edges.as_ref().unwrap();
        let vh_1 = &vh[ev[edge][1]];
        let vh_1_pos = vh_1.iter().position(|h| h.edge() == edge).expect("inconsistent vertex/edge topology");

        // Add new edge. This is unchecked, so we must fix the topology ourselves.
        let e_data = EdgeData { vertices: [new_v, ev[edge][1]], ..self.edge_data(edge) };
        let new_e = self.add_edge_unchecked(e_data, false);

        // Replace the half-edge
        let ev = self.edges_vertices.as_mut().unwrap();
        let vh = self.vertices_half_edges.as_mut().unwrap();
        vh[ev[edge][1]][vh_1_pos] = HalfEdge::new(new_e, true);
        
        // Connect the split edges together
        ev[edge][1] = new_v;
        vh[new_v] = vec![HalfEdge::new(edge, true), HalfEdge::new(new_e, false)];

        // If there are no faces, we are done.
        let fh = if let Some(fh) = self.faces_half_edges.as_mut() { fh } else { return (new_v, new_e) };
        let ec = self.edges_face_corners.as_mut().unwrap();
        // Unfortunately, if the mesh is oriented and the split edge is adjacent to only one face,
        // then the order of vh[new_v] matters.
        if ec[edge][0].len() == 1 && ec[edge][1].len() == 0 { vh[new_v].reverse() }
        
        // Update faces to use split edge
        let mut corners = ec[edge].iter().enumerate()
            .flat_map(|(i, cs)| cs.iter().map(move |&c| (c, i)))
            .collect::<Vec<_>>();
        // Go in reverse order so indices don't get messed up
        corners.sort();
        corners.reverse();
        let half_edges = new_e.split();
        for (c, i) in corners {
            let insert_index = c.corner() + (1 - i);
            // Fix face corner numbering
            for j in (insert_index..fh[c.face()].len()).rev() {
                let used_h = fh.at(FaceCorner(c.face(), j));
                let used_c = ec.at_mut(used_h).iter_mut().find(|c| c == &&FaceCorner(c.face(), j))
                    .expect("inconsisistent edge/face topology");
                used_c.1 += 1;
            }
            // And add the new half-edge into the face
            fh[c.face()].insert(insert_index, half_edges[i]);
            ec.at_mut(half_edges[i]).push(FaceCorner(c.face(), insert_index));
        }

        (new_v, new_e)
    }
}