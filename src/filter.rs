use std::mem;

use indexmap::{indexmap, IndexMap};
use nalgebra::{DMatrix, DVector, Scalar};
use serde_json::Value;
use typed_index_collections::TiVec;

use crate::fold::{CoordsRef, Edge, EdgeData, EdgesFaceCornersEx, Face, FaceCorner, FaceData, FacesCustom, FacesHalfEdges, Frame, HalfEdge, Vertex, VertexData};

#[derive(Clone, Copy)]
struct SwapRemove;
#[derive(Clone, Copy)]
struct ShiftRemove;
trait RemoveStrategy: Clone {
    fn index_map<K>(&self, vec_len: usize, key: K, k: K) -> Option<K> where
        K: PartialOrd + From<usize>,
        usize: From<K>;
    fn remove<K, T>(&self, vec: &mut TiVec<K, T>, key: K) -> T where usize: From<K>;
    fn remove_column<T: Scalar>(&self, vec: &mut DMatrix<T>, column: usize) -> DVector<T>;
}

impl RemoveStrategy for SwapRemove {
    fn index_map<K>(&self, vec_len: usize, key: K, k: K) -> Option<K> where
        K: PartialOrd + From<usize>,
        usize: From<K>
    {
        if k == key { None }
        else if k == K::from(vec_len - 1) { Some(key) }
        else { Some(k) }
    }

    fn remove<K, T>(&self, vec: &mut TiVec<K, T>, key: K) -> T where usize: From<K> { vec.swap_remove(key) }

    fn remove_column<T: Scalar>(&self, vec: &mut DMatrix<T>, column: usize) -> DVector<T> {
        let result = vec.column(column).into_owned();
        vec.swap_columns(column, vec.ncols() - 1);
        let ncols = vec.ncols();
        *vec = mem::take(vec).remove_column(ncols - 1);
        result
    }
}

impl RemoveStrategy for ShiftRemove {
    fn index_map<K>(&self, _vec_len: usize, key: K, k: K) -> Option<K> where
        K: PartialOrd + From<usize>,
        usize: From<K>
    {
        if k == key { None }
        else if k < key { Some(k) }
        else { Some(K::from(usize::from(k) - 1)) }
    }

    fn remove<K, T>(&self, vec: &mut TiVec<K, T>, key: K) -> T where usize: From<K> { vec.remove(key) }

    fn remove_column<T: Scalar>(&self, vec: &mut DMatrix<T>, column: usize) -> DVector<T> {
        let result = vec.column(column).into_owned();
        *vec = mem::take(vec).remove_column(column);
        result
    }
}

impl Frame {
    /// Gets a reference to the vertex coordinates, either exact or approximate.
    pub fn coords_ref(&self) -> Option<CoordsRef<'_>> {
        if let Some(coords) = self.vertices_coords_exact.as_ref() {
            Some(CoordsRef::Exact(coords))
        } else if let Some(coords) = self.vertices_coords_f64.as_ref() {
            Some(CoordsRef::Approx(coords))
        } else { None }
    }

    fn remap_vertices_refs(&mut self, new_indices: TiVec<Vertex, Option<Vertex>>) {
        if let Some(edges_vertices) = self.edges_vertices.as_mut() {
            for (e, vertices) in edges_vertices.iter_mut_enumerated() {
                *vertices = vertices
                    .map(|v| new_indices[v].unwrap_or_else(|| panic!("tried to remap vertex {v} to None while an edge {e} contains it")))
            }
        }
    }

    fn remap_edges_refs(&mut self, new_indices: TiVec<Edge, Option<Edge>>) {
        if let Some(faces_half_edges) = self.faces_half_edges.as_mut() {
            for (f, half_edges) in faces_half_edges.iter_mut_enumerated() {
                *half_edges = half_edges.drain(..)
                    .map(|h| new_indices[h.edge()].map(|e| HalfEdge::new(e, h.flip_bit()))
                        .unwrap_or_else(|| panic!("tried to remap half-edge {h} to None while a face {f} contains it")))
                    .collect::<Vec<_>>();
            }
        }

        if let Some(vertices_half_edges) = self.vertices_half_edges.as_mut() {
            for half_edges in vertices_half_edges {
                *half_edges = half_edges.drain(..)
                    .flat_map(|h| new_indices[h.edge()].map(|e| HalfEdge::new(e, h.flip_bit())))
                    .collect::<Vec<_>>();
            }
        }

        if let Some(edge_orders) = self.edge_orders.as_mut() {
            *edge_orders = edge_orders.drain(..)
                .flat_map(|(e1, e2, order)|
                    new_indices[e1].and_then(|e1| new_indices[e2].map(|e2| (e1, e2, order))))
                .collect::<Vec<_>>();
        }
    }

    fn remap_faces_refs(&mut self, new_indices: TiVec<Face, Option<Face>>) {
        if let Some(edges_face_corners) = self.edges_face_corners.as_mut() {
            for (_, corners) in edges_face_corners.half_iter_mut_enumerated() {
                *corners = corners.drain(..)
                    .flat_map(|c| new_indices[c.face()].map(|f| FaceCorner(f, c.corner())))
                    .collect::<Vec<_>>();
            }
        }

        if let Some(face_orders) = self.face_orders.as_mut() {
            *face_orders = face_orders.drain(..)
                .flat_map(|(f1, f2, order)|
                    new_indices[f1].and_then(|f1| new_indices[f2].map(|f2| (f1, f2, order))))
                .collect::<Vec<_>>();
        }
    }

    /// Removes all edges and faces attached to the vertex as well and returns their data.
    fn remove_vertex_generic(
        &mut self,
        vertex: Vertex,
        strategy: impl RemoveStrategy,
    ) -> (VertexData, IndexMap<Edge, EdgeData>, IndexMap<Face, FaceData>) {
        let face_data = if let Some(edges_face_corners) = self.edges_face_corners.as_ref() {
            let mut faces = self.vertices_half_edges.as_ref().unwrap()[vertex].iter()
                .flat_map(|&h| edges_face_corners.at(h).iter().map(|c| c.face()))
                .collect::<Vec<_>>();
            faces.sort();
            faces.dedup();
            faces.reverse();
            faces.into_iter().map(|f| (f, self.remove_face_generic(f, strategy.clone()))).collect::<IndexMap<_, _>>()
        } else { indexmap! {} };

        let edge_data = if let Some(vertices_half_edges) = self.vertices_half_edges.as_ref() {
            let mut edges = vertices_half_edges[vertex].iter().map(|h| h.edge()).collect::<Vec<_>>();
            edges.sort();
            edges.reverse(); // no need for dedup because self-loops aren't allowed
            edges.into_iter().map(|e| (e, self.remove_edge_generic(e, strategy.clone()).0)).collect::<IndexMap<_, _>>()
        } else { indexmap! {} };

        let new_indices = (0..self.num_vertices).map(
            |v| strategy.index_map(self.num_vertices, vertex, Vertex(v))).collect::<TiVec<_, _>>();
        self.num_vertices -= 1;
        
        let data = VertexData {
            coords_f64: self.vertices_coords_f64.as_mut().map(|vec| strategy.remove_column(vec, vertex.0)),
            coords_exact: self.vertices_coords_exact.as_mut().map(|vec| strategy.remove_column(vec, vertex.0)),
            half_edges: self.vertices_half_edges.as_mut().map(|vec| strategy.remove(vec, vertex)),
            custom: self.vertices_custom.iter_mut()
                .map(|(field, value)| (field.clone(), strategy.remove(value, vertex)))
                    .collect::<IndexMap<_, _>>()
        };
        self.remap_vertices_refs(new_indices);
        
        (data, edge_data, face_data)
    }

    /// Shift-removes a vertex `vertex`, shifting all vertexs after it to fill the gap,
    /// and shift-removes all edges and faces attached to the vertex.
    /// Returns the vertex data for that vertex, and the edge and face data for all removed edges and faces.
    pub fn shift_remove_vertex(&mut self, vertex: Vertex) -> (VertexData, IndexMap<Edge, EdgeData>, IndexMap<Face, FaceData>) {
        self.remove_vertex_generic(vertex, ShiftRemove)
    }

    /// Swap-removes a vertex `vertex`, moving the last vertex into its position,
    /// and swap-removes all edges and faces attached to the vertex.
    /// Returns the vertex data for that vertex, and the edge and face data for all removed edges and faces.
    pub fn swap_remove_vertex(&mut self, vertex: Vertex) -> (VertexData, IndexMap<Edge, EdgeData>, IndexMap<Face, FaceData>) {
        self.remove_vertex_generic(vertex, SwapRemove)
    }

    /// Removes all faces attached to the edge as well and returns their data.
    fn remove_edge_generic(
        &mut self,
        edge: Edge,
        strategy: impl RemoveStrategy,
    ) -> (EdgeData, IndexMap<Face, FaceData>) {
        let face_data = if let Some(edges_face_corners) = self.edges_face_corners.as_ref() {
            let mut faces = edges_face_corners[edge].iter().flatten().map(|c| c.face()).collect::<Vec<_>>();
            faces.sort();
            faces.dedup();
            faces.reverse();
            faces.into_iter().map(|f| (f, self.remove_face_generic(f, strategy.clone()))).collect::<IndexMap<_, _>>()
        } else { indexmap! {} };

        let edges_vertices = self.edges_vertices.as_mut().expect("tried to remove an edge without there being any");

        let new_indices = edges_vertices.iter_enumerated().map(
            |(e, _)| strategy.index_map(edges_vertices.len(), edge, e)).collect::<TiVec<_, _>>();
        
        let data = EdgeData {
            vertices: strategy.remove(edges_vertices, edge),
            face_corners: self.edges_face_corners.as_mut().map(|vec| strategy.remove(vec, edge)),
            assignment: self.edges_assignment.as_mut().map(|vec| strategy.remove(vec, edge)),
            fold_angle_f64: self.edges_fold_angle_f64.as_mut().map(|vec| strategy.remove(vec, edge)),
            fold_angle_exact: self.edges_fold_angle_exact.as_mut().map(|vec| strategy.remove(vec, edge)),
            length_f64: self.edges_length_f64.as_mut().map(|vec| strategy.remove(vec, edge)),
            length2_exact: self.edges_length2_exact.as_mut().map(|vec| strategy.remove(vec, edge)),
            custom: self.edges_custom.iter_mut()
                .map(|(field, value)| (field.clone(), strategy.remove(value, edge)))
                    .collect::<IndexMap<_, _>>()
        };
        self.remap_edges_refs(new_indices);
        
        (data, face_data)
    }

    /// Shift-removes a edge `edge`, shifting all edges after it to fill the gap,
    /// and shift-removes all faces attached to the edge.
    /// Returns the edge data for that edge, and the face data for all removed faces.
    pub fn shift_remove_edge(&mut self, edge: Edge) -> (EdgeData, IndexMap<Face, FaceData>) {
        self.remove_edge_generic(edge, ShiftRemove)
    }

    /// Swap-removes a edge `edge`, moving the last edge into its position,
    /// and swap-removes all faces attached to the edge.
    /// Returns the edge data for that edge, and the face data for all removed faces.
    pub fn swap_remove_edge(&mut self, edge: Edge) -> (EdgeData, IndexMap<Face, FaceData>) {
        self.remove_edge_generic(edge, SwapRemove)
    }

    fn remove_face_generic(
        &mut self,
        face: Face,
        strategy: impl RemoveStrategy,
    ) -> FaceData {
        let faces_half_edges = self.faces_half_edges.as_mut().expect("tried to remove a face without there being any");

        let new_indices = faces_half_edges.iter_enumerated().map(
            |(f, _)| strategy.index_map(faces_half_edges.len(), face, f)).collect::<TiVec<_, _>>();
        
        let half_edges = strategy.remove(faces_half_edges, face);
        let custom = self.faces_custom.iter_mut()
            .map(|(field, value)| (field.clone(), strategy.remove(value, face)))
            .collect::<IndexMap<_, _>>();
        self.remap_faces_refs(new_indices);
        
        FaceData { half_edges, custom }
    }

    /// Shift-removes a face `face`, shifting all faces after it to fill the gap.
    /// Returns the face data for that face.
    pub fn shift_remove_face(&mut self, face: Face) -> FaceData {
        self.remove_face_generic(face, ShiftRemove)
    }

    /// Swap-removes a face `face`, moving the last face into its position.
    /// Returns the face data for that face.
    pub fn swap_remove_face(&mut self, face: Face) -> FaceData {
        self.remove_face_generic(face, SwapRemove)
    }
}

/// Given two edges defined by their vertices, gets a vertex incident to both of them, if one exists.
pub fn edges_vertices_incident(e1: [Vertex; 2], e2: [Vertex; 2]) -> Option<Vertex> {
    e1.into_iter().flat_map(|v1| e2.into_iter().find(|v2| v1 == *v2)).next()
}

/// Given two vertices defined by their half-edges, gets the half-edges incident to both of them.
/// The half-edges point in the v1->v2 direction.
/// 
/// Make sure the vertices aren't the exact same vertex before calling this.
pub fn vertices_half_edges_incident(v1: &[HalfEdge], v2: &[HalfEdge]) -> Vec<HalfEdge> {
    v1.iter().copied().filter(|h1| v2.iter().find(|h2| h2.edge() == h1.edge()).is_some()).collect::<Vec<_>>()
}

/// Given an edge and a vertex on that edge, gets the other vertex
pub fn try_other_vertex(edge: [Vertex; 2], vertex: Vertex) -> Option<Vertex> {
    if edge[0] == vertex { Some(edge[1]) } else if edge[1] == vertex { Some(edge[0]) } else { None }
}

/// Given an edge and a vertex on that edge, gets the other vertex
pub fn other_vertex(edge: [Vertex; 2], vertex: Vertex) -> Vertex {
    if edge[0] == vertex { edge[1] } else { edge[0] }
}

#[cfg(test)]
mod test {
    use core::f64;

    use exact_number::based_expr;
    use indexmap::{IndexMap, indexmap};
    use nalgebra::{DMatrix, DVector};
    use serde_json::json;
    use typed_index_collections::ti_vec;

    use crate::{filter::{edges_vertices_incident, other_vertex}, fold::{EdgeAssignment, EdgeData, EdgeOrder, FaceData, FaceOrder, Frame, VertexData}, geom::ExactAngle};
    use crate::fold::{Vertex as V, Edge as E, Face as F, HalfEdge as H, FaceCorner as C};

    #[test]
    fn test_edges_vertices_incident() {
        assert_eq!(edges_vertices_incident([V(0), V(2)], [V(2), V(4)]), Some(V(2)));
        assert_eq!(edges_vertices_incident([V(0), V(2)], [V(4), V(2)]), Some(V(2)));
        assert_eq!(edges_vertices_incident([V(2), V(0)], [V(2), V(4)]), Some(V(2)));
        assert_eq!(edges_vertices_incident([V(2), V(0)], [V(4), V(2)]), Some(V(2)));
        assert_eq!(edges_vertices_incident([V(2), V(0)], [V(4), V(5)]), None);
    }

    #[test]
    fn test_other_vertex() {
        assert_eq!(other_vertex([V(5), V(4)], V(4)), V(5));
        assert_eq!(other_vertex([V(5), V(4)], V(5)), V(4));
        assert_eq!(other_vertex([V(3), V(7)], V(3)), V(7));
        assert_eq!(other_vertex([V(3), V(7)], V(7)), V(3));
    }

    #[test]
    fn test_swap_remove_vertex() {
        use EdgeAssignment::*;
        
        // 0---3---1---4---2
        // |       |       |
        // 0   0   1   1   2
        // |       |       |
        // 3---5---4---6---5
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            num_vertices: 6,
            vertices_coords_f64: Some(DMatrix::from_vec(2, 6, vec![
                0.0, 1.0,
                2.0, 1.0,
                4.0, 1.0,
                0.0, 0.0,
                2.0, 0.0,
                4.0, 0.0,
            ])),
            vertices_coords_exact: Some(DMatrix::from_vec(2, 6, vec![
                based_expr!(0), based_expr!(1),
                based_expr!(2), based_expr!(1),
                based_expr!(4), based_expr!(1),
                based_expr!(0), based_expr!(0),
                based_expr!(2), based_expr!(0),
                based_expr!(4), based_expr!(0),
            ])),
            vertices_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(3), json!(4), json!(5), json!(6)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(12), H(3), H(11)],
            vec![H(5), H(13)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
            vec![H(2), H(12), H(5), H(9)],
        ]));

        let mut frame = wrap_fold.clone();
        let (vertex_data, edge_datas, face_datas) = frame.swap_remove_vertex(V(0));
        assert_eq!(vertex_data, VertexData {
            coords_f64: Some(DVector::from_vec(vec![0.0, 1.0])),
            coords_exact: Some(DVector::from_vec(vec![based_expr!(0), based_expr!(1)])),
            half_edges: Some(vec![]),
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(edge_datas, indexmap! {
            E(0) => EdgeData {
                vertices: [V(0), V(3)],
                face_corners: Some([vec![], vec![]]),
                assignment: None,
                fold_angle_f64: None,
                fold_angle_exact: None,
                length_f64: None,
                length2_exact: None,
                custom: indexmap! {}
            },
            E(3) => EdgeData {
                vertices: [V(0), V(1)],
                face_corners: Some([vec![], vec![]]),
                assignment: None,
                fold_angle_f64: None,
                fold_angle_exact: None,
                length_f64: None,
                length2_exact: None,
                custom: indexmap! {}
            },
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.num_vertices, 5);
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(3), V(4)],
            [V(1), V(4)],
            [V(2), V(0)],
            [V(4), V(0)],
            [V(1), V(2)],
        ]));
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(5), H(7)],
            vec![H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(0)],
            vec![H(6), H(3), H(1)],
        ]));
        assert_eq!(frame.vertices_coords_f64, Some(DMatrix::from_vec(2, 5, vec![
            4.0, 0.0,
            2.0, 1.0,
            4.0, 1.0,
            0.0, 0.0,
            2.0, 0.0,
        ])));
        assert_eq!(frame.vertices_coords_exact, Some(DMatrix::from_vec(2, 5, vec![
            based_expr!(4), based_expr!(0),
            based_expr!(2), based_expr!(1),
            based_expr!(4), based_expr!(1),
            based_expr!(0), based_expr!(0),
            based_expr!(2), based_expr!(0),
        ])));
        assert_eq!(frame.vertices_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(6), json!(2), json!(3), json!(4), json!(5)]
        });
    }

    #[test]
    fn test_shift_remove_edge() {
        use EdgeAssignment::*;
        
        // 0---3---1---4---2
        // |       |       |
        // 0   0   1   1   2
        // |       |       |
        // 3---5---4---6---5
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            edges_assignment: Some(ti_vec![Boundary, Valley, Boundary, Boundary, Boundary, Boundary, Boundary]),
            edges_fold_angle_f64: Some(ti_vec![0.0, -180.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            edges_fold_angle_exact: Some(ti_vec![
                ExactAngle::ZERO, ExactAngle::NEG_PI, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO,
            ]),
            edges_length_f64: Some(ti_vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]),
            edges_length2_exact: Some(ti_vec![
                based_expr!(1), based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2), based_expr!(2)
            ]),
            edge_orders: Some(vec![
                (E(0), E(1), EdgeOrder::Right),
                (E(1), E(0), EdgeOrder::Left),
                (E(1), E(2), EdgeOrder::Right),
                (E(2), E(1), EdgeOrder::Left),
            ]),
            edges_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(3), json!(4), json!(5), json!(6), json!(7)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(12), H(3), H(11)],
            vec![H(5), H(13)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
            vec![H(2), H(12), H(5), H(9)],
        ]));

        let mut frame = wrap_fold.clone();
        let (edge_data, face_datas) = frame.shift_remove_edge(E(0));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(0), V(3)],
            face_corners: Some([vec![], vec![]]),
            assignment: Some(EdgeAssignment::Boundary),
            fold_angle_f64: Some(0.0),
            fold_angle_exact: Some(ExactAngle::ZERO),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(4)],
            vec![H(5), H(0), H(6)],
            vec![H(7), H(2)],
            vec![H(8)],
            vec![H(10), H(1), H(9)],
            vec![H(3), H(11)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(1), V(4)],
            [V(2), V(5)],
            [V(0), V(1)],
            [V(1), V(2)],
            [V(3), V(4)],
            [V(4), V(5)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 0)], vec![]],
            [vec![], vec![C(F(0), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![]],
            [vec![C(F(0), 1)], vec![]],
        ]));
        assert_eq!(frame.edges_assignment, Some(ti_vec![Valley, Boundary, Boundary, Boundary, Boundary, Boundary]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![-180.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![
            ExactAngle::NEG_PI, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO,
        ]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![1.0, 1.0, 2.0, 2.0, 2.0, 2.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![
            based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2), based_expr!(2)
        ]));
        assert_eq!(frame.edge_orders, Some(vec![
            (E(0), E(1), EdgeOrder::Right),
            (E(1), E(0), EdgeOrder::Left),
        ]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(2), json!(3), json!(4), json!(5), json!(6), json!(7)]
        });

        // Remove the edge adjacent to both faces
        let mut frame = wrap_fold.clone();
        let (edge_data, face_datas) = frame.shift_remove_edge(E(1));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(1), V(4)],
            face_corners: Some([vec![], vec![]]),
            assignment: Some(EdgeAssignment::Valley),
            fold_angle_f64: Some(-180.0),
            fold_angle_exact: Some(ExactAngle::NEG_PI),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(2) }
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            },
            F(1) => FaceData {
                half_edges: vec![H(2), H(12), H(5), H(9)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(4)],
            vec![H(5), H(6)],
            vec![H(7), H(2)],
            vec![H(8), H(1)],
            vec![H(10), H(9)],
            vec![H(3), H(11)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(3)],
            [V(2), V(5)],
            [V(0), V(1)],
            [V(1), V(2)],
            [V(3), V(4)],
            [V(4), V(5)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.edges_assignment, Some(ti_vec![Boundary, Boundary, Boundary, Boundary, Boundary, Boundary]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![
            ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO
        ]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![1.0, 1.0, 2.0, 2.0, 2.0, 2.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![
            based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2), based_expr!(2)
        ]));
        assert_eq!(frame.edge_orders, Some(vec![]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(3), json!(4), json!(5), json!(6), json!(7)]
        });

        // Slitted face, remove edge adjacent to the same face multiple times
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
        let (edge_data, face_datas) = frame.shift_remove_edge(E(8));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(0), V(4)],
            face_corners: Some([vec![], vec![]]),
            assignment: None,
            fold_angle_f64: None,
            fold_angle_exact: None,
            length_f64: None,
            length2_exact: None,
            custom: indexmap! {}
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
                custom: indexmap! {}
            },
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(8), H(10), H(12), H(14)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(4), V(5)],
            [V(5), V(6)],
            [V(6), V(7)],
            [V(7), V(4)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(0), 2)], vec![]],
            [vec![C(F(0), 3)], vec![]],
        ]));
    }

    #[test]
    fn test_swap_remove_edge() {
        use EdgeAssignment::*;
        
        // 0---3---1---4---2
        // |       |       |
        // 0   0   1   1   2
        // |       |       |
        // 3---5---4---6---5
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            edges_assignment: Some(ti_vec![Boundary, Valley, Boundary, Boundary, Boundary, Boundary, Boundary]),
            edges_fold_angle_f64: Some(ti_vec![0.0, -180.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            edges_fold_angle_exact: Some(ti_vec![
                ExactAngle::ZERO, ExactAngle::NEG_PI, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO,
            ]),
            edges_length_f64: Some(ti_vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]),
            edges_length2_exact: Some(ti_vec![
                based_expr!(1), based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2), based_expr!(2)
            ]),
            edge_orders: Some(vec![
                (E(0), E(1), EdgeOrder::Right),
                (E(1), E(0), EdgeOrder::Left),
                (E(1), E(2), EdgeOrder::Right),
                (E(2), E(1), EdgeOrder::Left),
            ]),
            edges_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(3), json!(4), json!(5), json!(6), json!(7)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(12), H(3), H(11)],
            vec![H(5), H(13)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
            vec![H(2), H(12), H(5), H(9)],
        ]));

        let mut frame = wrap_fold.clone();
        let (edge_data, face_datas) = frame.swap_remove_edge(E(0));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(0), V(3)],
            face_corners: Some([vec![], vec![]]),
            assignment: Some(EdgeAssignment::Boundary),
            fold_angle_f64: Some(0.0),
            fold_angle_exact: Some(ExactAngle::ZERO),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10)],
            vec![H(0), H(3), H(11)],
            vec![H(5), H(1)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(2), H(0), H(5), H(9)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(4), V(5)],
            [V(1), V(4)],
            [V(2), V(5)],
            [V(0), V(1)],
            [V(1), V(2)],
            [V(3), V(4)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(0), 0)], vec![]],
            [vec![], vec![C(F(0), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.edges_assignment, Some(ti_vec![Boundary, Valley, Boundary, Boundary, Boundary, Boundary]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![0.0, -180.0, 0.0, 0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![
            ExactAngle::ZERO, ExactAngle::NEG_PI, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO
        ]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![2.0, 1.0, 1.0, 2.0, 2.0, 2.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![
            based_expr!(2), based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2)
        ]));
        assert_eq!(frame.edge_orders, Some(vec![
            (E(1), E(2), EdgeOrder::Right),
            (E(2), E(1), EdgeOrder::Left),
        ]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(7), json!(2), json!(3), json!(4), json!(5), json!(6)]
        });

        // Remove the edge adjacent to both faces
        let mut frame = wrap_fold.clone();
        let (edge_data, face_datas) = frame.swap_remove_edge(E(1));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(1), V(4)],
            face_corners: Some([vec![], vec![]]),
            assignment: Some(EdgeAssignment::Valley),
            fold_angle_f64: Some(-180.0),
            fold_angle_exact: Some(ExactAngle::NEG_PI),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(2) }
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            },
            F(1) => FaceData {
                half_edges: vec![H(2), H(12), H(5), H(9)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(2),  H(11)],
            vec![H(5), H(3)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(3)],
            [V(4), V(5)],
            [V(2), V(5)],
            [V(0), V(1)],
            [V(1), V(2)],
            [V(3), V(4)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.edges_assignment, Some(ti_vec![Boundary, Boundary, Boundary, Boundary, Boundary, Boundary]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![
            ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO, ExactAngle::ZERO
        ]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![1.0, 2.0, 1.0, 2.0, 2.0, 2.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![
            based_expr!(1), based_expr!(2), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2)
        ]));
        assert_eq!(frame.edge_orders, Some(vec![]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(7), json!(3), json!(4), json!(5), json!(6)]
        });

        // Slitted face, remove edge adjacent to the same face multiple times
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
        let (edge_data, face_datas) = frame.swap_remove_edge(E(8));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(0), V(4)],
            face_corners: Some([vec![], vec![]]),
            assignment: None,
            fold_angle_f64: None,
            fold_angle_exact: None,
            length_f64: None,
            length2_exact: None,
            custom: indexmap! {}
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
                custom: indexmap! {}
            },
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(8), H(10), H(12), H(14)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(4), V(5)],
            [V(5), V(6)],
            [V(6), V(7)],
            [V(7), V(4)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(0), 2)], vec![]],
            [vec![C(F(0), 3)], vec![]],
        ]));
    }

    #[test]
    fn test_shift_remove_face() {
        // 0---4---1---5---2---6---3
        // |       |       |       |
        // 0   0   1   1   2   2   3
        // |       |       |       |
        // 4---7---5---8---6---9---7
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            face_orders: Some(vec![
                (F(0), F(1), FaceOrder::Above),
                (F(1), F(0), FaceOrder::Above),
                (F(0), F(2), FaceOrder::Above),
                (F(2), F(0), FaceOrder::Below),
                (F(2), F(1), FaceOrder::Above),
                (F(1), F(2), FaceOrder::Above),
            ]),
            faces_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(4)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(8)],
            vec![H(9), H(2), H(10)],
            vec![H(11), H(4), H(12)],
            vec![H(13), H(6)],
            vec![H(14), H(1)],
            vec![H(16), H(3), H(15)],
            vec![H(18), H(5), H(17)],
            vec![H(7), H(19)],
        ]), Some(ti_vec![
            vec![H(0), H(14), H(3), H(9)],
            vec![H(2), H(16), H(5), H(11)],
            vec![H(4), H(18), H(7), H(13)],
        ]));

        let mut frame = wrap_fold.clone();
        let face_data = frame.shift_remove_face(F(0));
        assert_eq!(face_data, FaceData {
            half_edges: vec![H(0), H(14), H(3), H(9)],
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(2), H(16), H(5), H(11)],
            vec![H(4), H(18), H(7), H(13)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(1), 0)], vec![C(F(0), 2)]],
            [vec![], vec![C(F(1), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![C(F(1), 3)]],
            [vec![], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(1), 1)], vec![]],
        ]));
        assert_eq!(frame.face_orders, Some(vec![
            (F(1), F(0), FaceOrder::Above),
            (F(0), F(1), FaceOrder::Above),
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(2), json!(4)]
        });

        // Removing the last entry
        let mut frame = wrap_fold.clone();
        let face_data = frame.shift_remove_face(F(2));
        assert_eq!(face_data, FaceData {
            half_edges: vec![H(4), H(18), H(7), H(13)],
            custom: indexmap! { "test:test".to_owned() => json!(4) }
        });
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(14), H(3), H(9)],
            vec![H(2), H(16), H(5), H(11)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(1), 0)], vec![C(F(0), 2)]],
            [vec![], vec![C(F(1), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![C(F(1), 3)]],
            [vec![], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(1), 1)], vec![]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.face_orders, Some(vec![
            (F(0), F(1), FaceOrder::Above),
            (F(1), F(0), FaceOrder::Above),
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(2)]
        });
    }

    #[test]
    fn test_swap_remove_face() {
        // 0---4---1---5---2---6---3
        // |       |       |       |
        // 0   0   1   1   2   2   3
        // |       |       |       |
        // 4---7---5---8---6---9---7
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            face_orders: Some(vec![
                (F(0), F(1), FaceOrder::Above),
                (F(1), F(0), FaceOrder::Above),
                (F(0), F(2), FaceOrder::Above),
                (F(2), F(0), FaceOrder::Below),
                (F(2), F(1), FaceOrder::Above),
                (F(1), F(2), FaceOrder::Above),
            ]),
            faces_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(4)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(8)],
            vec![H(9), H(2), H(10)],
            vec![H(11), H(4), H(12)],
            vec![H(13), H(6)],
            vec![H(14), H(1)],
            vec![H(16), H(3), H(15)],
            vec![H(18), H(5), H(17)],
            vec![H(7), H(19)],
        ]), Some(ti_vec![
            vec![H(0), H(14), H(3), H(9)],
            vec![H(2), H(16), H(5), H(11)],
            vec![H(4), H(18), H(7), H(13)],
        ]));

        let mut frame = wrap_fold.clone();
        let face_data = frame.swap_remove_face(F(0));
        assert_eq!(face_data, FaceData {
            half_edges: vec![H(0), H(14), H(3), H(9)],
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(4), H(18), H(7), H(13)],
            vec![H(2), H(16), H(5), H(11)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![C(F(1), 0)], vec![]],
            [vec![C(F(0), 0)], vec![C(F(1), 2)]],
            [vec![], vec![C(F(0), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(1), 3)]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![]],
            [vec![C(F(1), 1)], vec![]],
            [vec![C(F(0), 1)], vec![]],
        ]));
        assert_eq!(frame.face_orders, Some(vec![
            (F(0), F(1), FaceOrder::Above),
            (F(1), F(0), FaceOrder::Above),
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(4), json!(2)]
        });

        // Removing the last entry
        let mut frame = wrap_fold.clone();
        let face_data = frame.swap_remove_face(F(2));
        assert_eq!(face_data, FaceData {
            half_edges: vec![H(4), H(18), H(7), H(13)],
            custom: indexmap! { "test:test".to_owned() => json!(4) }
        });
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(14), H(3), H(9)],
            vec![H(2), H(16), H(5), H(11)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(1), 0)], vec![C(F(0), 2)]],
            [vec![], vec![C(F(1), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![C(F(1), 3)]],
            [vec![], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(1), 1)], vec![]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.face_orders, Some(vec![
            (F(0), F(1), FaceOrder::Above),
            (F(1), F(0), FaceOrder::Above),
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(2)]
        });
    }
}