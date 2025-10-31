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

#[cfg(test)]
mod test {
    use indexmap::indexset;
    use nalgebra::{DVector, vector};
    use typed_index_collections::ti_vec;

    use crate::{Frame, FrameAttribute};
    use crate::{Vertex as V, Edge as E, HalfEdge as H, FaceCorner as C, Face as F};

    #[test]
    fn test_split_edge() {
        // This is annoying to test; especially because you can turn connected components inside-out,
        // so just print results.
 
        // No faces. Just a 3-star.
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ]), None);
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(0), DVector::from_vec(vec![]));
        assert_eq!(new_vertex, V(4));
        assert_eq!(new_edge, E(3));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(0)], [V(0), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(1)]);
        frame.assert_topology_consistent();

        // A simple triangle.
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ]), Some(ti_vec![
            vec![H(5), H(3), H(1)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(1), DVector::from_vec(vec![]));
        assert_eq!(new_vertex, V(3));
        assert_eq!(new_edge, E(3));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(1)], [V(1), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(2)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(5), new_edge.split()[1], H(3), H(1)]);
        frame.assert_topology_consistent();

        // Now it's oriented the other way. The half-edge order better be right.
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(5)],
            vec![H(2), H(1)],
            vec![H(4), H(3)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(1), DVector::from_vec(vec![]));
        assert_eq!(new_vertex, V(3));
        assert_eq!(new_edge, E(3));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(1)], [V(1), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(2)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), H(2), new_edge.split()[0], H(4)]);
        frame.assert_topology_consistent();

        // A doubled square
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)]
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(3), DVector::from_vec(vec![]));
        let [hfwd, hbck] = new_edge.split();
        assert_eq!(new_vertex, V(4));
        assert_eq!(new_edge, E(4));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(3)], [V(3), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(0)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), H(2), H(4), H(6), hfwd]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(1)], vec![hbck, H(7), H(5), H(3), H(1)]);
        frame.assert_topology_consistent();

        // A slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let mut frame = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(16)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15), H(17)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
            vec![H(8), H(10), H(12), H(14)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(8), DVector::from_vec(vec![]));
        let [hfwd, hbck] = new_edge.split();
        assert_eq!(new_vertex, V(8));
        assert_eq!(new_edge, E(9));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(8)], [V(0), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(4)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), H(2), H(4), H(6), H(16), hfwd, H(15), H(13), H(11), H(9), hbck, H(17)]);
        frame.assert_topology_consistent();

        // A möbius face. Noticable, a face duplicates a half-edge.
        // Edges go 0->1->2->3->0 and 0->2.
        // »---0       3--3»
        //     |"0   ."     
        //     4  "."       
        //     | 2" ".      
        // >---2"     "1--1>  
        let mut frame = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9), H(7), H(5), H(9)],
            vec![H(0), H(2), H(4), H(6)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(4), DVector::from_vec(vec![]));
        let [hfwd, hbck] = new_edge.split();
        assert_eq!(new_vertex, V(4));
        assert_eq!(new_edge, E(5));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(4)], [V(0), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(2)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), H(2), hbck, H(9), H(7), H(5), hbck, H(9)]);
        frame.assert_topology_consistent();

        // Three triangles from a single edge
        let mut frame = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(5), H(9), H(13)],
            vec![H(1), H(2), H(6), H(10)],
            vec![H(3), H(4)],
            vec![H(7), H(8)],
            vec![H(11), H(12)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(0), H(6), H(8)],
            vec![H(0), H(10), H(12)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(0), DVector::from_vec(vec![]));
        let [hfwd, hbck] = new_edge.split();
        assert_eq!(new_vertex, V(5));
        assert_eq!(new_edge, E(7));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(0)], [V(0), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(1)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), hfwd, H(2), H(4)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(1)], vec![H(0), hfwd, H(6), H(8)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(2)], vec![H(0), hfwd, H(10), H(12)]);
        frame.assert_topology_consistent();
    }
}