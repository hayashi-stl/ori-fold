use indexmap::IndexMap;
use nalgebra::{ClosedSubAssign, DMatrix, Scalar, Vector2};

use crate::{filter, fold::{CoordsRef, EdgeAssignment, Fold, Frame}, geom::{sort_by_angle_ref, AngleRep}};
use crate::geom;

/// Given `edges_vertices` (defining edge endpoints),
/// automatically computes the `vertices_edges` property.
/// However, note that the `vertices_edges` arrays will not be sorted in counterclockwise order.
pub fn edges_vertices_to_vertices_edges_unsorted(edges_vertices: &[[usize; 2]], num_vertices: usize) -> Vec<Vec<usize>> {
    let mut vertices_edges = (0..num_vertices).map(|_| vec![]).collect::<Vec<_>>();
    for (edge, vertices) in edges_vertices.iter().copied().enumerate() {
        for vertex in vertices {
            vertices_edges[vertex].push(edge);
        }
    }
    vertices_edges
}

#[must_use]
pub fn sort_vertices_edges_generic<T>(mut vertices_edges: Vec<Vec<usize>>, edges_vertices: &[[usize; 2]], coords: &DMatrix<T>) -> Vec<Vec<usize>> where 
    T: Scalar + ClosedSubAssign,
    Vector2<T>: AngleRep,
{
    for (i, edges) in vertices_edges.iter_mut().enumerate() {
        geom::sort_by_angle_ref(edges, &i,
            |e| coords.fixed_view::<2, 1>(0, filter::other_vertex(edges_vertices[*e], i)));
    }
    vertices_edges
}

/// Given `vertices_edges` and `edges_vertices`
/// and some 2D `vertices_coords` (exact or approximate),
/// sorts each vertices_edges array in counterclockwise order around the vertex in the plane,
/// according to `AngleRep::angle_rep`, and returns the result.
#[must_use]
pub fn sort_vertices_edges(vertices_edges: Vec<Vec<usize>>, edges_vertices: &[[usize; 2]], coords: CoordsRef) -> Vec<Vec<usize>> {
    match coords {
        CoordsRef::Exact(c) => sort_vertices_edges_generic(vertices_edges, edges_vertices, c),
        CoordsRef::Approx(c) => sort_vertices_edges_generic(vertices_edges, edges_vertices, c),
    }
}

/// Given some 2D `vertices_coords` (exact or approximate)
/// and `edges_vertices` (defining edge endpoints),
/// automatically computes the `vertices_edges` properties and
/// sorts them counterclockwise by angle in the plane.
/// according to `AngleRep::angle_rep()`.
pub fn edges_vertices_to_vertices_edges_sorted(coords: CoordsRef, edges_vertices: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let vertices_edges = edges_vertices_to_vertices_edges_unsorted(edges_vertices, coords.num_vertices());
    sort_vertices_edges(vertices_edges, edges_vertices, coords)
}

/// Given `faces_edges` and `edges_vertices`, computes
/// `edges_faces`.
/// This assumes the frame is an orientable manifold,
/// and stores each entry of `edges_faces` as [left face, right face] when viewed from the edge,
/// with `None` used for faces that don't exist (boundary edges)
pub fn faces_edges_to_edges_faces_orientable(faces_edges: &[Vec<Vec<usize>>], edges_vertices: &[[usize; 2]]) -> Vec<Vec<Option<usize>>> {
    let edges_faces = vec![vec![None; 2]; edges_vertices.len()];
}

impl Frame {
    /// Given a FOLD object frame with `edges_vertices == Some(_)` (defining edge endpoints),
    /// automatically computes the `vertices_edges` property.
    /// However, note that the `vertices_edges` arrays will not be sorted in counterclockwise order.
    pub fn edges_vertices_to_vertices_edges_unsorted(&mut self) {
        let num_vertices = self.num_vertices;
        self.vertices_edges = Some(edges_vertices_to_vertices_edges_unsorted(
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
            num_vertices
        ));
    }

    /// Given a FOLD object with `vertices_coords_f64` or `vertices_coords_exact` is `Some(2D vectors)`
    /// and `vertices_edges == Some(_)`,
    /// sorts each vertices_edges array in counterclockwise order around the vertex in the plane,
    /// according to `AngleRep::angle_rep`.
    pub fn sort_vertices_edges(&mut self) {
        self.vertices_edges = Some(sort_vertices_edges(
            self.vertices_edges.take().expect("vertices_edges must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
            self.coords_ref().expect("vertices_coords_exact or vertices_coords_f64 must exist")
        ));
    }

    /// Given a FOLD object frame with `vertices_coords_f64` or `vertices_coords_exact` is `Some(2D vectors)`
    /// and `edges_vertices == Some(_)` (defining edge endpoints),
    /// automatically computes the `vertices_vertices` and `vertices_edges` properties and
    /// sorts them counterclockwise by angle in the plane.
    /// according to `AngleRep::angle_rep()`.
    pub fn edges_vertices_to_vertices_edges_sorted(&mut self) {
        self.vertices_edges = Some(edges_vertices_to_vertices_edges_sorted(
            self.coords_ref().expect("vertices_coords_exact or vertices_coords_f64 must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
        ));
    }
}

#[cfg(test)]
mod test {
    use exact_number::based_expr;
    use nalgebra::DMatrix;

    use crate::fold::Frame;

    #[test]
    fn test_edges_vertices_to_vertices_edges_no_edges() {
        let mut frame = Frame {
            edges_vertices: Some(vec![]),
            ..Default::default()
        };
        frame.edges_vertices_to_vertices_edges_unsorted();
        assert_eq!(frame.vertices_edges, Some(vec![]));
    }

    #[test]
    fn test_edges_vertices_to_vertices_edges_star() {
        let mut frame = Frame {
            edges_vertices: Some(vec![
                [0, 1],
                [0, 2],
                [0, 3],
            ]),
            ..Default::default()
        };
        frame.edges_vertices_to_vertices_edges_unsorted();
        frame.vertices_edges.as_mut().unwrap()[0].sort();
        assert_eq!(frame.vertices_edges, Some(vec![
            vec![0, 1, 2],
            vec![0],
            vec![1],
            vec![2],
        ]));
    }

    #[test]
    fn test_edges_vertices_to_vertices_edges_square_slash() {
        let mut frame = Frame {
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [1, 3],
            ]),
            ..Default::default()
        };
        frame.edges_vertices_to_vertices_edges_unsorted();
        frame.vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(frame.vertices_edges, Some(vec![
            vec![0, 3],
            vec![0, 1, 4],
            vec![1, 2],
            vec![2, 3, 4],
        ]));
    }

    #[test]
    fn test_edges_vertices_to_vertices_edges_isolated_vertex() {
        let mut frame = Frame {
            vertices_coords_f64: Some(DMatrix::from_vec(2, 4, vec![
                0.0, 1.0,
                2.0, 3.0,
                4.0, 1.0,
                6.0, 8.0,
            ])),
            edges_vertices: Some(vec![
                [0, 1],
                [0, 2],
            ]),
            ..Default::default()
        };
        frame.edges_vertices_to_vertices_edges_unsorted();
        frame.vertices_edges.as_mut().unwrap()[0].sort();
        assert_eq!(frame.vertices_edges, Some(vec![
            vec![0, 1],
            vec![0],
            vec![1],
            vec![],
        ]));
    }

    #[test]
    fn test_sort_vertices_vertices_nothing() {
        let mut frame = Frame {
            vertices_coords_f64: Some(DMatrix::repeat(2, 0, 0.0)),
            vertices_vertices: Some(vec![]),
            ..Default::default()
        };
        frame.sort_vertices_vertices();
        assert_eq!(frame.vertices_vertices, Some(vec![]));
    }

    #[test]
    fn test_sort_vertices_vertices_triangle() {
        let mut frame = Frame {
            vertices_coords_f64: Some(DMatrix::from_vec(2, 3, vec![
                0.0, 0.0,
                2.0, 0.0,
                1.0, 1.0
            ])),
            vertices_vertices: Some(vec![
                vec![1, 2],
                vec![0, 2],
                vec![0, 1]
            ]),
            ..Default::default()
        };
        frame.sort_vertices_vertices();
        assert_eq!(frame.vertices_vertices, Some(vec![
            vec![1, 2],
            vec![2, 0],
            vec![0, 1]
        ]));
    }

    #[test]
    fn test_sort_vertices_vertices_triangle_prefer_exact() {
        let mut frame = Frame {
            vertices_coords_f64: Some(DMatrix::from_vec(2, 3, vec![
                0.0, 0.0,
                -2.0, 0.0,
                -1.0, -1.0 // Intentional wrong coordinates
            ])),
            vertices_coords_exact: Some(DMatrix::from_vec(2, 3, vec![
                based_expr!(0), based_expr!(0),
                based_expr!(2), based_expr!(0),
                based_expr!(1), based_expr!(1),
            ])),
            vertices_vertices: Some(vec![
                vec![1, 2],
                vec![0, 2],
                vec![0, 1]
            ]),
            ..Default::default()
        };
        frame.sort_vertices_vertices();
        assert_eq!(frame.vertices_vertices, Some(vec![
            vec![1, 2],
            vec![2, 0],
            vec![0, 1]
        ]));
    }

    #[test]
    fn test_vertices_vertices_to_vertices_edges_single_edge() {
        let mut frame = Frame {
            vertices_vertices: Some(vec![
                vec![1],
                vec![0],
            ]),
            edges_vertices: Some(vec![
                [0, 1],
            ]),
            ..Default::default()
        };
        frame.vertices_vertices_to_vertices_edges();
        assert_eq!(frame.vertices_edges, Some(vec![
            vec![0],
            vec![0],
        ]))
    }

    #[test]
    fn test_vertices_vertices_to_vertices_edges_square_cross() {
        let mut frame = Frame {
            vertices_vertices: Some(vec![
                vec![1, 4, 3],
                vec![2, 4, 0],
                vec![3, 4, 1],
                vec![0, 4, 2],
                vec![0, 1, 2, 3]
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [0, 4],
                [1, 4],
                [2, 4],
                [3, 4]
            ]),
            ..Default::default()
        };
        frame.vertices_vertices_to_vertices_edges();
        assert_eq!(frame.vertices_edges, Some(vec![
            vec![0, 4, 3],
            vec![1, 5, 0],
            vec![2, 6, 1],
            vec![3, 7, 2],
            vec![4, 5, 6, 7]
        ]))
    }

    #[test]
    fn test_vertices_vertices_to_faces_vertices_nothing() {
        let mut frame = Frame {
            vertices_vertices: Some(vec![]),
            ..Default::default()
        };
        frame.vertices_vertices_to_faces_vertices();
        assert_eq!(frame.faces_vertices, Some(vec![]));
    }

    #[test]
    fn test_vertices_vertices_to_faces_vertices_single_face() {
        let mut frame = Frame {
            vertices_vertices: Some(vec![
                vec![1, 2],
                vec![2, 0],
                vec![0, 1]
            ]),
            ..Default::default()
        };
        frame.vertices_vertices_to_faces_vertices();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 0) { faces_vertices[1].rotate_left(k) }
        faces_vertices.sort();
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 2],
            vec![0, 2, 1]
        ]);
    }

    #[test]
    fn test_vertices_vertices_to_faces_vertices_slitted_face() {
        // 3---------2
        // |         |
        // |  7---6  |
        // |  |   |  |
        // |  4---5  |
        // |,"       |
        // 0---------1
        let mut frame = Frame {
            vertices_vertices: Some(vec![
                vec![1, 4, 3],
                vec![2, 0],
                vec![3, 1],
                vec![0, 2],
                vec![5, 7, 0],
                vec![6, 4],
                vec![7, 5],
                vec![4, 6]
            ]),
            ..Default::default()
        };
        frame.vertices_vertices_to_faces_vertices();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        if let Some((i, _)) = faces_vertices[0].iter().enumerate().max_by_key(|(_, k)| **k) { faces_vertices[0].rotate_left(i) }
        if let Some((i, _)) = faces_vertices[1].iter().enumerate().max_by_key(|(_, k)| **k) { faces_vertices[1].rotate_left(i) }
        if let Some((i, _)) = faces_vertices[2].iter().enumerate().max_by_key(|(_, k)| **k) { faces_vertices[2].rotate_left(i) }
        faces_vertices.sort();
        assert_eq!(faces_vertices, &vec![
            vec![3, 2, 1, 0],
            vec![7, 4, 5, 6],
            vec![7, 6, 5, 4, 0, 1, 2, 3, 0, 4],
        ]);
    }

    #[test]
    fn test_vertices_vertices_to_faces_vertices_square_cross() {
        let mut frame = Frame {
            vertices_vertices: Some(vec![
                vec![1, 4, 3],
                vec![2, 4, 0],
                vec![3, 4, 1],
                vec![0, 4, 2],
                vec![0, 1, 2, 3]
            ]),
            ..Default::default()
        };
        frame.vertices_vertices_to_faces_vertices();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        faces_vertices.sort_by_key(|vs| (vs.len(), vs.contains(&0) && vs.contains(&3), vs.iter().sum::<usize>()));
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 1) { faces_vertices[1].rotate_left(k) }
        if let Some(k) = faces_vertices[2].iter().position(|x| *x == 2) { faces_vertices[2].rotate_left(k) }
        if let Some(k) = faces_vertices[3].iter().position(|x| *x == 3) { faces_vertices[3].rotate_left(k) }
        if let Some(k) = faces_vertices[4].iter().position(|x| *x == 3) { faces_vertices[4].rotate_left(k) }
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
            vec![3, 2, 1, 0],
        ]);
    }

    #[test]
    fn test_vertices_vertices_to_faces_vertices_square_cross_3d_coords() {
        let mut frame = Frame {
            vertices_coords_exact: Some(DMatrix::from_vec(3, 5, vec![
                based_expr!(0), based_expr!(0), based_expr!(0),
                based_expr!(1), based_expr!(0), based_expr!(0),
                based_expr!(1), based_expr!(1), based_expr!(0),
                based_expr!(0), based_expr!(1), based_expr!(0),
                based_expr!(1/2), based_expr!(1/2), based_expr!(0),
            ])),
            vertices_vertices: Some(vec![
                vec![1, 4, 3],
                vec![2, 4, 0],
                vec![3, 4, 1],
                vec![0, 4, 2],
                vec![0, 1, 2, 3]
            ]),
            ..Default::default()
        };
        frame.vertices_vertices_to_faces_vertices();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        faces_vertices.sort_by_key(|vs| (vs.len(), vs.contains(&0) && vs.contains(&3), vs.iter().sum::<usize>()));
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 1) { faces_vertices[1].rotate_left(k) }
        if let Some(k) = faces_vertices[2].iter().position(|x| *x == 2) { faces_vertices[2].rotate_left(k) }
        if let Some(k) = faces_vertices[3].iter().position(|x| *x == 3) { faces_vertices[3].rotate_left(k) }
        if let Some(k) = faces_vertices[4].iter().position(|x| *x == 3) { faces_vertices[4].rotate_left(k) }
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
            vec![3, 2, 1, 0],
        ]);
    }

    #[test]
    fn test_vertices_vertices_to_faces_vertices_square_cross_2d_coords() {
        let mut frame = Frame {
            vertices_coords_exact: Some(DMatrix::from_vec(2, 5, vec![
                based_expr!(0), based_expr!(0),
                based_expr!(1), based_expr!(0),
                based_expr!(1), based_expr!(1),
                based_expr!(0), based_expr!(1),
                based_expr!(1/2), based_expr!(1/2),
            ])),
            vertices_vertices: Some(vec![
                vec![1, 4, 3],
                vec![2, 4, 0],
                vec![3, 4, 1],
                vec![0, 4, 2],
                vec![0, 1, 2, 3]
            ]),
            ..Default::default()
        };
        frame.vertices_vertices_to_faces_vertices();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        faces_vertices.sort_by_key(|vs| (vs.len(), vs.contains(&0) && vs.contains(&3), vs.iter().sum::<usize>()));
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 1) { faces_vertices[1].rotate_left(k) }
        if let Some(k) = faces_vertices[2].iter().position(|x| *x == 2) { faces_vertices[2].rotate_left(k) }
        if let Some(k) = faces_vertices[3].iter().position(|x| *x == 3) { faces_vertices[3].rotate_left(k) }
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
        ]);
    }

    #[test]
    fn test_vertices_edges_to_faces_vertices_edges_multi_edge() {
        let mut frame = Frame {
            vertices_edges: Some(vec![
                vec![0, 2, 3],
                vec![1, 0],
                vec![3, 2, 1]
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [0, 2],
                [0, 2],
            ]),
            ..Default::default()
        };
        frame.vertices_edges_to_faces_vertices_edges();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        let faces_edges = frame.faces_edges.as_mut().unwrap();
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 0) { faces_vertices[1].rotate_left(k) }
        if let Some(k) = faces_vertices[2].iter().position(|x| *x == 0) { faces_vertices[2].rotate_left(k) }
        faces_vertices.sort();
        faces_edges.sort_by_key(|es| es.len());
        if let Some(k) = faces_edges[0].iter().position(|x| *x == 2) { faces_edges[0].rotate_left(k) }
        if let Some(k) = faces_edges[1].iter().position(|x| *x == 0) { faces_edges[1].rotate_left(k) }
        if let Some(k) = faces_edges[2].iter().position(|x| *x == 0) { faces_edges[2].rotate_left(k) }
        faces_edges.sort();
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 2],
            vec![0, 2],
            vec![0, 2, 1],
        ]);
        assert_eq!(faces_edges, &vec![
            vec![0, 1, 2],
            vec![0, 3, 1],
            vec![2, 3],
        ]);
    }

    #[test]
    fn test_vertices_edges_to_faces_vertices_edges_single_face() {
        let mut frame = Frame {
            vertices_edges: Some(vec![
                vec![0, 2],
                vec![1, 0],
                vec![2, 1]
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [0, 2],
            ]),
            ..Default::default()
        };
        frame.vertices_edges_to_faces_vertices_edges();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        let faces_edges = frame.faces_edges.as_mut().unwrap();
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 0) { faces_vertices[1].rotate_left(k) }
        faces_vertices.sort();
        if let Some(k) = faces_edges[0].iter().position(|x| *x == 0) { faces_edges[0].rotate_left(k) }
        if let Some(k) = faces_edges[1].iter().position(|x| *x == 0) { faces_edges[1].rotate_left(k) }
        faces_edges.sort();
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 2],
            vec![0, 2, 1],
        ]);
        assert_eq!(faces_edges, &vec![
            vec![0, 1, 2],
            vec![0, 2, 1],
        ]);
    }

    #[test]
    fn test_vertices_edges_to_faces_vertices_edges_slitted_face() {
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let mut frame = Frame {
            vertices_edges: Some(vec![
                vec![0, 8, 3],
                vec![1, 0],
                vec![2, 1],
                vec![3, 2],
                vec![4, 7, 8],
                vec![5, 4],
                vec![6, 5],
                vec![7, 6],
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [4, 5],
                [5, 6],
                [6, 7],
                [7, 4],
                [0, 4],
            ]),
            ..Default::default()
        };
        frame.vertices_edges_to_faces_vertices_edges();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        let faces_edges = frame.faces_edges.as_mut().unwrap();
        if let Some((i, _)) = faces_vertices[0].iter().enumerate().max_by_key(|(_, k)| **k) { faces_vertices[0].rotate_left(i) }
        if let Some((i, _)) = faces_vertices[1].iter().enumerate().max_by_key(|(_, k)| **k) { faces_vertices[1].rotate_left(i) }
        if let Some((i, _)) = faces_vertices[2].iter().enumerate().max_by_key(|(_, k)| **k) { faces_vertices[2].rotate_left(i) }
        faces_vertices.sort();
        if let Some((i, _)) = faces_edges[0].iter().enumerate().min_by_key(|(_, k)| **k) { faces_edges[0].rotate_left(i) }
        if let Some((i, _)) = faces_edges[1].iter().enumerate().min_by_key(|(_, k)| **k) { faces_edges[1].rotate_left(i) }
        if let Some((i, _)) = faces_edges[2].iter().enumerate().min_by_key(|(_, k)| **k) { faces_edges[2].rotate_left(i) }
        faces_edges.sort();
        assert_eq!(faces_vertices, &vec![
            vec![3, 2, 1, 0],
            vec![7, 4, 5, 6],
            vec![7, 6, 5, 4, 0, 1, 2, 3, 0, 4],
        ]);
        assert_eq!(faces_edges, &vec![
            vec![0, 1, 2, 3, 8, 7, 6, 5, 4, 8],
            vec![0, 3, 2, 1],
            vec![4, 5, 6, 7],
        ]);
    }

    #[test]
    fn test_vertices_edges_to_faces_vertices_edges_square_cross() {
        let mut frame = Frame {
            vertices_edges: Some(vec![
                vec![0, 4, 3],
                vec![1, 5, 0],
                vec![2, 6, 1],
                vec![3, 7, 2],
                vec![4, 5, 6, 7]
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [0, 4],
                [1, 4],
                [2, 4],
                [3, 4],
            ]),
            ..Default::default()
        };
        frame.vertices_edges_to_faces_vertices_edges();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        let faces_edges = frame.faces_edges.as_mut().unwrap();
        faces_vertices.sort_by_key(|vs| (vs.len(), vs.contains(&0) && vs.contains(&3), vs.iter().sum::<usize>()));
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 1) { faces_vertices[1].rotate_left(k) }
        if let Some(k) = faces_vertices[2].iter().position(|x| *x == 2) { faces_vertices[2].rotate_left(k) }
        if let Some(k) = faces_vertices[3].iter().position(|x| *x == 3) { faces_vertices[3].rotate_left(k) }
        if let Some(k) = faces_vertices[4].iter().position(|x| *x == 3) { faces_vertices[4].rotate_left(k) }
        faces_edges.sort_by_key(|vs| (vs.len(), vs.contains(&4) && vs.contains(&7), vs.iter().sum::<usize>()));
        if let Some(k) = faces_edges[0].iter().position(|x| *x == 0) { faces_edges[0].rotate_left(k) }
        if let Some(k) = faces_edges[1].iter().position(|x| *x == 1) { faces_edges[1].rotate_left(k) }
        if let Some(k) = faces_edges[2].iter().position(|x| *x == 2) { faces_edges[2].rotate_left(k) }
        if let Some(k) = faces_edges[3].iter().position(|x| *x == 3) { faces_edges[3].rotate_left(k) }
        if let Some(k) = faces_edges[4].iter().position(|x| *x == 3) { faces_edges[4].rotate_left(k) }
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
            vec![3, 2, 1, 0],
        ]);
        assert_eq!(faces_edges, &vec![
            vec![0, 5, 4],
            vec![1, 6, 5],
            vec![2, 7, 6],
            vec![3, 4, 7],
            vec![3, 2, 1, 0],
        ]);
    }

    #[test]
    fn test_vertices_edges_to_faces_vertices_edges_square_cross_3d_coords() {
        let mut frame = Frame {
            vertices_coords_exact: Some(DMatrix::from_vec(3, 5, vec![
                based_expr!(0), based_expr!(0), based_expr!(0),
                based_expr!(1), based_expr!(0), based_expr!(0),
                based_expr!(1), based_expr!(1), based_expr!(0),
                based_expr!(0), based_expr!(1), based_expr!(0),
                based_expr!(1/2), based_expr!(1/2), based_expr!(0),
            ])),
            vertices_edges: Some(vec![
                vec![0, 4, 3],
                vec![1, 5, 0],
                vec![2, 6, 1],
                vec![3, 7, 2],
                vec![4, 5, 6, 7]
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [0, 4],
                [1, 4],
                [2, 4],
                [3, 4],
            ]),
            ..Default::default()
        };
        frame.vertices_edges_to_faces_vertices_edges();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        let faces_edges = frame.faces_edges.as_mut().unwrap();
        faces_vertices.sort_by_key(|vs| (vs.len(), vs.contains(&0) && vs.contains(&3), vs.iter().sum::<usize>()));
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 1) { faces_vertices[1].rotate_left(k) }
        if let Some(k) = faces_vertices[2].iter().position(|x| *x == 2) { faces_vertices[2].rotate_left(k) }
        if let Some(k) = faces_vertices[3].iter().position(|x| *x == 3) { faces_vertices[3].rotate_left(k) }
        if let Some(k) = faces_vertices[4].iter().position(|x| *x == 3) { faces_vertices[4].rotate_left(k) }
        faces_edges.sort_by_key(|vs| (vs.len(), vs.contains(&4) && vs.contains(&7), vs.iter().sum::<usize>()));
        if let Some(k) = faces_edges[0].iter().position(|x| *x == 0) { faces_edges[0].rotate_left(k) }
        if let Some(k) = faces_edges[1].iter().position(|x| *x == 1) { faces_edges[1].rotate_left(k) }
        if let Some(k) = faces_edges[2].iter().position(|x| *x == 2) { faces_edges[2].rotate_left(k) }
        if let Some(k) = faces_edges[3].iter().position(|x| *x == 3) { faces_edges[3].rotate_left(k) }
        if let Some(k) = faces_edges[4].iter().position(|x| *x == 3) { faces_edges[4].rotate_left(k) }
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
            vec![3, 2, 1, 0],
        ]);
        assert_eq!(faces_edges, &vec![
            vec![0, 5, 4],
            vec![1, 6, 5],
            vec![2, 7, 6],
            vec![3, 4, 7],
            vec![3, 2, 1, 0],
        ]);
    }

    #[test]
    fn test_vertices_edges_to_faces_vertices_edges_square_cross_2d_coords() {
        let mut frame = Frame {
            vertices_coords_exact: Some(DMatrix::from_vec(2, 5, vec![
                based_expr!(0), based_expr!(0),
                based_expr!(1), based_expr!(0),
                based_expr!(1), based_expr!(1),
                based_expr!(0), based_expr!(1),
                based_expr!(1/2), based_expr!(1/2),
            ])),
            vertices_edges: Some(vec![
                vec![0, 4, 3],
                vec![1, 5, 0],
                vec![2, 6, 1],
                vec![3, 7, 2],
                vec![4, 5, 6, 7]
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [0, 4],
                [1, 4],
                [2, 4],
                [3, 4],
            ]),
            ..Default::default()
        };
        frame.vertices_edges_to_faces_vertices_edges();
        let faces_vertices = frame.faces_vertices.as_mut().unwrap();
        let faces_edges = frame.faces_edges.as_mut().unwrap();
        faces_vertices.sort_by_key(|vs| (vs.len(), vs.contains(&0) && vs.contains(&3), vs.iter().sum::<usize>()));
        if let Some(k) = faces_vertices[0].iter().position(|x| *x == 0) { faces_vertices[0].rotate_left(k) }
        if let Some(k) = faces_vertices[1].iter().position(|x| *x == 1) { faces_vertices[1].rotate_left(k) }
        if let Some(k) = faces_vertices[2].iter().position(|x| *x == 2) { faces_vertices[2].rotate_left(k) }
        if let Some(k) = faces_vertices[3].iter().position(|x| *x == 3) { faces_vertices[3].rotate_left(k) }
        faces_edges.sort_by_key(|vs| (vs.len(), vs.contains(&4) && vs.contains(&7), vs.iter().sum::<usize>()));
        if let Some(k) = faces_edges[0].iter().position(|x| *x == 0) { faces_edges[0].rotate_left(k) }
        if let Some(k) = faces_edges[1].iter().position(|x| *x == 1) { faces_edges[1].rotate_left(k) }
        if let Some(k) = faces_edges[2].iter().position(|x| *x == 2) { faces_edges[2].rotate_left(k) }
        if let Some(k) = faces_edges[3].iter().position(|x| *x == 3) { faces_edges[3].rotate_left(k) }
        assert_eq!(faces_vertices, &vec![
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
        ]);
        assert_eq!(faces_edges, &vec![
            vec![0, 5, 4],
            vec![1, 6, 5],
            vec![2, 7, 6],
            vec![3, 4, 7],
        ]);
    }

    #[test]
    fn test_edges_vertices_to_edges_faces_edges_square_cross() {
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 4],
                vec![1, 2, 4],
                vec![2, 3, 4],
                vec![3, 0, 4]
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [0, 4],
                [1, 4],
                [2, 4],
                [3, 4]
            ]),
            ..Default::default()
        };
        frame.edges_vertices_to_edges_faces_edges();
        assert_eq!(frame.edges_faces, Some(vec![
            vec![Some(0), None],
            vec![Some(1), None],
            vec![Some(2), None],
            vec![Some(3), None],
            vec![Some(3), Some(0)],
            vec![Some(0), Some(1)],
            vec![Some(1), Some(2)],
            vec![Some(2), Some(3)]
        ]));
        assert_eq!(frame.faces_edges, Some(vec![
            vec![0, 5, 4],
            vec![1, 6, 5],
            vec![2, 7, 6],
            vec![3, 4, 7]
        ]));
    }

    #[test]
    fn test_edges_vertices_to_edges_faces_edges_single_face() {
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 2]
            ]),
            edges_vertices: Some(vec![
                [1, 2],
                [0, 2],
                [0, 1],
            ]),
            ..Default::default()
        };
        frame.edges_vertices_to_edges_faces_edges();
        assert_eq!(frame.edges_faces, Some(vec![
            vec![Some(0), None],
            vec![None, Some(0)],
            vec![Some(0), None],
        ]));
        assert_eq!(frame.faces_edges, Some(vec![
            vec![2, 0, 1],
        ]));
    }

    #[test]
    fn test_edges_vertices_to_edges_faces_edges_slitted_face() {
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 2, 3, 0, 4, 7, 6, 5, 4],
                vec![4, 5, 6, 7],
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [4, 5],
                [5, 6],
                [6, 7],
                [7, 4],
                [0, 4],
            ]),
            ..Default::default()
        };
        frame.edges_vertices_to_edges_faces_edges();
        assert_eq!(frame.edges_faces, Some(vec![
            vec![Some(0), None],
            vec![Some(0), None],
            vec![Some(0), None],
            vec![Some(0), None],
            vec![Some(1), Some(0)],
            vec![Some(1), Some(0)],
            vec![Some(1), Some(0)],
            vec![Some(1), Some(0)],
            vec![Some(0), Some(0)],
        ]));
        assert_eq!(frame.faces_edges, Some(vec![
            vec![0, 1, 2, 3, 8, 7, 6, 5, 4, 8],
            vec![4, 5, 6, 7],
        ]));
    }

    #[test]
    fn test_faces_vertices_to_faces_edges_square_cross() {
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 4],
                vec![1, 2, 4],
                vec![2, 3, 4],
                vec![3, 0, 4]
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [0, 4],
                [1, 4],
                [2, 4],
                [3, 4]
            ]),
            ..Default::default()
        };
        frame.faces_vertices_to_faces_edges();
        assert_eq!(frame.faces_edges, Some(vec![
            vec![0, 5, 4],
            vec![1, 6, 5],
            vec![2, 7, 6],
            vec![3, 4, 7]
        ]));
    }

    #[test]
    fn test_faces_vertices_to_faces_edges_single_face() {
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 2]
            ]),
            edges_vertices: Some(vec![
                [1, 2],
                [0, 2],
                [0, 1],
            ]),
            ..Default::default()
        };
        frame.faces_vertices_to_faces_edges();
        assert_eq!(frame.faces_edges, Some(vec![
            vec![2, 0, 1],
        ]));
    }

    #[test]
    fn test_faces_vertices_to_faces_edges_slitted_face() {
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 2, 3, 0, 4, 7, 6, 5, 4],
                vec![4, 5, 6, 7],
            ]),
            edges_vertices: Some(vec![
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [4, 5],
                [5, 6],
                [6, 7],
                [7, 4],
                [0, 4],
            ]),
            ..Default::default()
        };
        frame.faces_vertices_to_faces_edges();
        assert_eq!(frame.faces_edges, Some(vec![
            vec![0, 1, 2, 3, 8, 7, 6, 5, 4, 8],
            vec![4, 5, 6, 7],
        ]));
    }

    #[test]
    fn test_faces_vertices_to_edges_square_cross() {
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 4],
                vec![1, 2, 4],
                vec![2, 3, 4],
                vec![3, 0, 4]
            ]),
            ..Default::default()
        };
        frame.faces_vertices_to_edges();
        // Make sure the result is reasonable; too lazy to canonicalize edges
        println!("faces_vertices_to_edges square cross:");
        println!("    edges_vertices: {:?}", frame.edges_vertices);
        println!("    edges_faces: {:?}", frame.edges_faces);
        println!("    faces_edges: {:?}", frame.faces_edges);
        println!("    edges_assignment: {:?}", frame.edges_assignment);
    }

    #[test]
    fn test_faces_vertices_to_edges_single_face() {
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 2]
            ]),
            ..Default::default()
        };
        frame.faces_vertices_to_edges();
        // Make sure the result is reasonable; too lazy to canonicalize edges
        println!("faces_vertices_to_edges single face:");
        println!("    edges_vertices: {:?}", frame.edges_vertices);
        println!("    edges_faces: {:?}", frame.edges_faces);
        println!("    faces_edges: {:?}", frame.faces_edges);
        println!("    edges_assignment: {:?}", frame.edges_assignment);
    }

    #[test]
    fn test_faces_vertices_to_edges_slitted_face() {
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let mut frame = Frame {
            faces_vertices: Some(vec![
                vec![0, 1, 2, 3, 0, 4, 7, 6, 5, 4],
                vec![4, 5, 6, 7],
            ]),
            ..Default::default()
        };
        frame.faces_vertices_to_edges();
        // Make sure the result is reasonable; too lazy to canonicalize edges
        println!("faces_vertices_to_edges slitted face:");
        println!("    edges_vertices: {:?}", frame.edges_vertices);
        println!("    edges_faces: {:?}", frame.edges_faces);
        println!("    faces_edges: {:?}", frame.faces_edges);
        println!("    edges_assignment: {:?}", frame.edges_assignment);
    }
}