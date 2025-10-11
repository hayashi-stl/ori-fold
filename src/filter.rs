use crate::fold::{CoordsRef, Frame};

impl Frame {
    /// Counts the number of vertices
    /// Looks in fields labelled `vertices_*` (reserved and custom),
    /// then  in fields labelled `*_vertices` (reserved and custom)
    /// 
    /// Returns 0 if it doesn't find a vertex field.
    pub fn num_vertices(&self) -> usize {
        vec![
            self.vertices_coords_f64.as_ref().map(|v| v.ncols()),
            self.vertices_coords_exact.as_ref().map(|v| v.ncols()),
            self.vertices_vertices.as_ref().map(|v| v.len()),
            self.vertices_edges.as_ref().map(|v| v.len()),
            self.vertices_faces.as_ref().map(|v| v.len()),
        ].into_iter().flatten().next()
            .or_else(|| self.vertices_custom.iter().map(|(_, v)| v.len()).next())
            .or_else(|| vec![
                self.edges_vertices.as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| n + 1),
                self.faces_vertices.as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| n + 1),
            ].into_iter().flatten().max())
            .unwrap_or(0)
    }

    /// Gets a reference to the vertex coordinates, either exact or approximate.
    pub fn coords_ref(&self) -> Option<CoordsRef<'_>> {
        if let Some(coords) = self.vertices_coords_exact.as_ref() {
            Some(CoordsRef::Exact(coords))
        } else if let Some(coords) = self.vertices_coords_f64.as_ref() {
            Some(CoordsRef::Approx(coords))
        } else { None }
    }
}

/// Given two edges defined by their vertices, gets a vertex incident to both of them, if one exists.
pub fn edges_vertices_incident(e1: [usize; 2], e2: [usize; 2]) -> Option<usize> {
    e1.into_iter().flat_map(|v1| e2.into_iter().find(|v2| v1 == *v2)).next()
}

/// Given an edge and a vertex on that edge, gets the other vertex
pub fn other_vertex(edge: [usize; 2], vertex: usize) -> usize {
    if edge[0] == vertex { edge[1] } else { edge[0] }
}

#[cfg(test)]
mod test {
    use exact_number::based_expr;
    use indexmap::IndexMap;
    use nalgebra::DMatrix;
    use serde_json::json;

    use crate::{filter::{edges_vertices_incident, other_vertex}, fold::Frame};

    #[test]
    fn test_num_vertices_vertices_coords_f64() {
        let frame = Frame {
            vertices_coords_f64: Some(DMatrix::from_vec(2, 2, vec![
                1.0, 2.0,
                3.0, 4.0
            ])),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 2);
    }

    #[test]
    fn test_num_vertices_vertices_coords_exact() {
        let frame = Frame {
            vertices_coords_exact: Some(DMatrix::from_vec(2, 1, vec![
                based_expr!(1), based_expr!(2),
            ])),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 1);
    }

    #[test]
    fn test_num_vertices_vertices_vertices() {
        let frame = Frame {
            vertices_vertices: Some(vec![
                vec![1, 4, 5],
                vec![3, 4, 1, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 2);
    }

    #[test]
    fn test_num_vertices_vertices_edges() {
        let frame = Frame {
            vertices_edges: Some(vec![
                vec![1, 4, 5],
                vec![3, 4, 1, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 2);
    }

    #[test]
    fn test_num_vertices_vertices_faces() {
        let frame = Frame {
            vertices_faces: Some(vec![
                vec![Some(1), None, Some(3)],
                vec![Some(0), Some(3), Some(4)],
                vec![None, Some(2), Some(5)],
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 3);
    }

    #[test]
    fn test_num_vertices_vertices_custom() {
        let frame = Frame {
            vertices_custom: vec![
                ("test:field".to_owned(), vec![
                    json!(3),
                    json!(5),
                    json!(5),
                ])
            ].into_iter().collect::<IndexMap<_, _>>(),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 3);
    }

    #[test]
    fn test_num_vertices_edges_vertices() {
        let frame = Frame {
            edges_vertices: Some(vec![
                [3, 4],
                [1, 4],
                [5, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 6);
    }

    #[test]
    fn test_num_vertices_faces_vertices() {
        let frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 3, 4],
                vec![2, 1, 4],
                vec![3, 5, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 6);
    }

    #[test]
    fn test_num_vertices_edges_vertices_faces_vertices() {
        let mut frame = Frame {
            edges_vertices: Some(vec![
                [1, 4],
                [2, 3],
            ]),
            faces_vertices: Some(vec![
                vec![2, 7, 4],
                vec![3, 4, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 8);
        frame.edges_vertices = Some(vec![
            [10, 2],
            [3, 5]
        ]);
        assert_eq!(frame.num_vertices(), 11);
    }

    #[test]
    fn test_num_vertices_no_vertex_field() {
        let mut frame = Frame {
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(), 0);
    }

    #[test]
    fn test_edges_vertices_incident() {
        assert_eq!(edges_vertices_incident([0, 2], [2, 4]), Some(2));
        assert_eq!(edges_vertices_incident([0, 2], [4, 2]), Some(2));
        assert_eq!(edges_vertices_incident([2, 0], [2, 4]), Some(2));
        assert_eq!(edges_vertices_incident([2, 0], [4, 2]), Some(2));
        assert_eq!(edges_vertices_incident([2, 0], [4, 5]), None);
    }

    #[test]
    fn test_other_vertex() {
        assert_eq!(other_vertex([5, 4], 4), 5);
        assert_eq!(other_vertex([5, 4], 5), 4);
        assert_eq!(other_vertex([3, 7], 3), 7);
        assert_eq!(other_vertex([3, 7], 7), 3);
    }
}