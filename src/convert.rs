use indexmap::IndexMap;
use nalgebra::{ClosedSubAssign, DMatrix, Scalar, Vector2};

use crate::{filter, fold::{CoordsRef, EdgeAssignment, Fold, Frame}, geom::{sort_by_angle_ref, AngleRep}};
use crate::geom;

/// Given `edges_vertices` (defining edge endpoints),
/// automatically computes the `vertices_vertices` property.
/// However, note that the `vertices_vertices` arrays will not be sorted in counterclockwise order.
pub fn edges_vertices_to_vertices_vertices_unsorted(edges_vertices: &[[usize; 2]], num_vertices: usize) -> Vec<Vec<usize>> {
    let mut vertices_vertices = (0..num_vertices).map(|_| vec![]).collect::<Vec<_>>();
    for [v, w] in edges_vertices.iter().copied() {
        vertices_vertices[v].push(w);
        vertices_vertices[w].push(v);
    }
    vertices_vertices
}

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
pub fn sort_vertices_vertices_generic<T>(mut vertices_vertices: Vec<Vec<usize>>, coords: &DMatrix<T>) -> Vec<Vec<usize>> where 
    T: Scalar + ClosedSubAssign,
    Vector2<T>: AngleRep,
{
    for (i, neighbors) in vertices_vertices.iter_mut().enumerate() {
        geom::sort_by_angle_ref(neighbors, &i, |x| coords.fixed_view::<2, 1>(0, *x));
    }
    vertices_vertices
}

/// Given `vertices_vertices`
/// and some 2D `vertices_coords` (exact or approximate),
/// sorts each vertices_vertices array in counterclockwise order around the vertex in the plane,
/// according to `AngleRep::angle_rep`, and returns the result.
#[must_use]
pub fn sort_vertices_vertices(vertices_vertices: Vec<Vec<usize>>, coords: CoordsRef) -> Vec<Vec<usize>> {
    match coords {
        CoordsRef::Exact(c) => sort_vertices_vertices_generic(vertices_vertices, c),
        CoordsRef::Approx(c) => sort_vertices_vertices_generic(vertices_vertices, c),
    }
}

/// Given some 2D `vertices_coords` (exact or approximate)
/// and `edges_vertices` (defining edge endpoints),
/// automatically computes the `vertices_vertices` property and
/// sorts them counterclockwise by angle in the plane.
/// according to `AngleRep::angle_rep()`.
pub fn edges_vertices_to_vertices_vertices_sorted(coords: CoordsRef, edges_vertices: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let result = edges_vertices_to_vertices_vertices_unsorted(edges_vertices, coords.num_vertices());
    sort_vertices_vertices(result, coords)
}

/// Given some 2D `vertices_coords` (exact or approximate)
/// and `edges_vertices` (defining edge endpoints),
/// automatically computes the `vertices_vertices` and `vertices_edges` properties and
/// sorts them counterclockwise by angle in the plane.
/// according to `AngleRep::angle_rep()`.
pub fn edges_vertices_to_vertices_edges_sorted(coords: CoordsRef, edges_vertices: &[[usize; 2]]) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let vertices_vertices = edges_vertices_to_vertices_vertices_sorted(coords, edges_vertices);
    let vertices_edges = vertices_vertices_to_vertices_edges(&vertices_vertices, edges_vertices);
    (vertices_vertices, vertices_edges)
}

/// Given `vertices_vertices` and `edges_vertices`,
/// calculates the corresponding `vertices_edges` property (preserving order).
pub fn vertices_vertices_to_vertices_edges(vertices_vertices: &[Vec<usize>], edges_vertices: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let edge_map = edges_vertices.iter().enumerate()
        .flat_map(|(i, [v1, v2])| [
            ((*v1, *v2), i),
            ((*v2, *v1), i)
        ]).collect::<IndexMap<_, _>>();
    let vertices_edges = vertices_vertices.iter().enumerate()
        .map(|(vertex, vertices)| vertices.iter().map(|v| edge_map[&(vertex, *v)]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    vertices_edges
}

/// Given counterclockwise-sorted `vertices_vertices`,
/// constructs the implicitly defined faces, setting the `faces_vertices` property.
/// If `vertices_coords` is `Some(2D points)`, then
/// this excludes faces that are not counterclockwise.
pub fn vertices_vertices_to_faces_vertices(vertices_vertices: &[Vec<usize>], coords: Option<CoordsRef>) -> Vec<Vec<usize>> {
    let mut next = vertices_vertices.iter().enumerate()
        .flat_map(|(v, neighbors)| neighbors.iter().enumerate()
            .map(move |(i, u)| ((*u, v), neighbors[(i + neighbors.len() - 1) % neighbors.len()])))
        .collect::<IndexMap<_, _>>();
    let mut faces_vertices = vec![];

    //for ((u, v), w) in &next {
    //    println!("{u}->{v}->{w}");
    //}
    while let Some(((mut u, mut v), mut w)) = next.pop() {
        let mut face = vec![u, v];
        while w != face[0] {
            face.push(w);
            (u, v) = (v, w);
            w = if let Some(w) = next.swap_remove(&(u, v)) { w } else {
                eprintln!("Confusion with face {face:?}");
                break;
            };
        }

        next.swap_remove(&(*face.last().unwrap(), face[0]));
        // Outside face is clockwise; exclude it. (if 2D and coordinates are defined)
        let exclude = match coords {
            None => false,
            Some(CoordsRef::Exact(c)) => c.nrows() == 2 &&
                    geom::polygon_orientation(DMatrix::from_columns(&face.iter().map(|v| c.column(*v)).collect::<Vec<_>>()).as_view()) <= 0,
            Some(CoordsRef::Approx(c)) => c.nrows() == 2 &&
                    geom::polygon_orientation(DMatrix::from_columns(&face.iter().map(|v| c.column(*v)).collect::<Vec<_>>()).as_view()) <= 0,
        };
        if !exclude {
            faces_vertices.push(face);
        }
    }
    faces_vertices
}

/// Given counterclockwise-sorted `vertices_edges` property,
/// and `edges_vertices`,
/// constructs the implicitly defined faces, setting both `faces_vertices` and `faces_edges` properties.
/// Handles multiple edges between the same two vertices (unlike `convert::vertices_vertices_to_faces_vertices`).
/// If `vertices_coords` is `Some(2D points)`, then
/// this excludes faces that are not counterclockwise.
pub fn vertices_edges_to_faces_vertices_edges(vertices_edges: &[Vec<usize>], edges_vertices: &[[usize; 2]], coords: Option<CoordsRef>)
    -> (Vec<Vec<usize>>, Vec<Vec<usize>>)
{
    let mut next = vertices_edges.iter().enumerate()
        .flat_map(|(v, neighbors)| neighbors.iter().enumerate()
            .map(|(i, e)| ((*e, v), neighbors[(i + neighbors.len() - 1) % neighbors.len()]))
            .collect::<Vec<_>>())
        .collect::<IndexMap<_, _>>();

    //for ((e1, v), e2) in &next {
    //    println!("{e1}->{v}->{e2}");
    //}
    let mut faces_vertices = vec![];
    let mut faces_edges = vec![];
    while let Some(((e1, vertex), mut e2)) = next.pop() {
        let mut edges = vec![e1];
        let mut vertices = vec![vertex];
        //println!("start edges: {edges:?}, vertices: {vertices:?}");
        while e2 != edges[0] {
            edges.push(e2);
            vertices.push(filter::other_vertex(edges_vertices[e2], *vertices.last().unwrap()));
            //println!("edges: {edges:?}, vertices: {vertices:?}");
            e2 = if let Some(e) = next.swap_remove(&(e2, *vertices.last().unwrap())) { e } else {
                eprintln!("Confusion with face containing edges {edges:?}");
                break;
            };
        }

        // Move e1 to the end so that edges[0] connects vertices[0] to vertices[1]
        edges.rotate_left(1);
        let exclude = match coords {
            None => false,
            Some(CoordsRef::Exact(c)) => c.nrows() == 2 &&
                    geom::polygon_orientation(DMatrix::from_columns(&vertices.iter().map(|v| c.column(*v)).collect::<Vec<_>>()).as_view()) <= 0,
            Some(CoordsRef::Approx(c)) => c.nrows() == 2 &&
                    geom::polygon_orientation(DMatrix::from_columns(&vertices.iter().map(|v| c.column(*v)).collect::<Vec<_>>()).as_view()) <= 0,
        };
        if !exclude {
            faces_vertices.push(vertices);
            faces_edges.push(edges);
        }
    }
    (faces_vertices, faces_edges)
}

/// Given some 2D `vertices_coords` (exact or approximate)
/// and `edges_vertices`,
/// computes a counterclockwise-sorted `vertices_vertices` property
/// and constructs the implicitly defined faces, setting `faces_vertices` property.
pub fn edges_vertices_to_faces_vertices(coords: CoordsRef, edges_vertices: &[[usize; 2]]) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let vertices_vertices = edges_vertices_to_vertices_vertices_sorted(coords, edges_vertices);
    let faces_vertices = vertices_vertices_to_faces_vertices(&vertices_vertices, Some(coords));
    (vertices_vertices, faces_vertices)
}

/// Given some 2D `vertices_coords` (exact or approximate)
/// and `edges_vertices`,
/// computes counterclockwise-sorted `vertices_vertices` and `vertices_edges` properties and
/// constructs the implicitly defined faces, setting both `faces_vertices` and `faces_edges` properties.
pub fn edges_vertices_to_faces_vertices_edges(coords: CoordsRef, edges_vertices: &[[usize; 2]])
    -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>)
{
    let (vertices_vertices, vertices_edges) = edges_vertices_to_vertices_edges_sorted(coords, edges_vertices);
    let (faces_vertices, faces_edges) = vertices_edges_to_faces_vertices_edges(&vertices_edges, edges_vertices, Some(coords));
    (vertices_vertices, vertices_edges, faces_vertices, faces_edges)
}

/// Given `faces_vertices`
/// and `edges_vertices`,
/// fills in the corresponding faces_edges property (preserving order).
pub fn faces_vertices_to_faces_edges(faces_vertices: &[Vec<usize>], edges_vertices: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let edge_map = edges_vertices.iter().enumerate()
        .flat_map(|(i, [v1, v2])| [
            ((*v1, *v2), i),
            ((*v2, *v1), i)
        ]).collect::<IndexMap<_, _>>();
    let faces_edges = faces_vertices.iter()
        .map(|vertices|
            vertices.iter().zip(vertices.iter().cycle().skip(1)).map(|(v0, v1)| edge_map[&(*v0, *v1)]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    faces_edges
}

/// Given `faces_vertices`, automatically computes
/// `edges_vertices`, `edges_faces`, `faces_edges`, and `edges_assignment`.
/// (indicating which edges are boundary with `EdgeAssignment::Boundary` and marking
/// all other edges with `EdgeAssignment::Unassigned`).
/// Currently assumes an orientable manifold, and uses `None`s to
/// represent missing neighbor faces in `edges_faces` (for boundary edges).
pub fn faces_vertices_to_edges(faces_vertices: &[Vec<usize>])
    -> (Vec<[usize; 2]>, Vec<Vec<Option<usize>>>, Vec<Vec<usize>>, Vec<EdgeAssignment>)
{
    let mut edges_vertices = vec![];
    let mut edges_faces = vec![];
    let mut faces_edges = vec![];
    let mut edges_assignment = vec![];
    let mut edge_map = IndexMap::<[usize; 2], usize>::new();
    
    for (face, vertices) in faces_vertices.iter().enumerate() {
        faces_edges.push({
            vertices.iter().enumerate().map(|(i, &v1)| {
                let v2 = vertices[(i + 1) % vertices.len()];
                let mut key = [v1, v2];
                key.sort();
                
                let edge = if let Some(&edge) = edge_map.get(&key) {
                    edges_assignment[edge] = EdgeAssignment::Unassigned;
                    edge
                } else {
                    let edge = edges_vertices.len();
                    edge_map.insert(key, edge);
                    edges_vertices.push(key);
                    edges_faces.push(vec![None; 2]);
                    edges_assignment.push(EdgeAssignment::Boundary);
                    edge
                };
                edges_faces[edge][if v1 <= v2 { 0 } else { 1 }] = Some(face);
                edge
            }).collect::<Vec<_>>()
        });
    }

    (edges_vertices, edges_faces, faces_edges, edges_assignment)
}

/// Given a FOLD object frame with `edges_vertices` and `faces_vertices`,
/// fills in `edges_faces` and `faces_edges`.
/// Currently assumes an orientable manifold, and uses `None`s to
/// represent missing neighbor faces in `edges_faces` (for boundary edges).
pub fn edges_vertices_to_edges_faces_edges(edges_vertices: &[[usize; 2]], faces_vertices: &[Vec<usize>])
    -> (Vec<Vec<Option<usize>>>, Vec<Vec<usize>>)
{
    let mut edges_faces = vec![vec![None, None]; edges_vertices.len()];
    let edge_map = edges_vertices.iter().enumerate().flat_map(|(edge, [v0, v1])| {[
        ((*v0, *v1), (edge, 0)),
        ((*v1, *v0), (edge, 1)),
    ]}).collect::<IndexMap<_, _>>();
    let faces_edges = faces_vertices.iter().enumerate().map(|(face, vertices)| {
        vertices.iter().enumerate().map(|(i, &v1)| {
            let v2 = vertices[(i + 1) % vertices.len()];
            let (edge, orient) = edge_map[&(v1, v2)];
            edges_faces[edge][orient] = Some(face);
            edge
        }).collect::<Vec<_>>()
    }).collect::<Vec<_>>();
    (edges_faces, faces_edges)
}

impl Frame {
    /// Given a FOLD object frame with `edges_vertices == Some(_)` (defining edge endpoints),
    /// automatically computes the `vertices_vertices` property.
    /// However, note that the `vertices_vertices` arrays will not be sorted in counterclockwise order.
    pub fn edges_vertices_to_vertices_vertices_unsorted(&mut self) {
        let num_vertices = self.num_vertices();
        self.vertices_vertices = Some(edges_vertices_to_vertices_vertices_unsorted(
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
            num_vertices
        ));
    }

    /// Given a FOLD object frame with `edges_vertices == Some(_)` (defining edge endpoints),
    /// automatically computes the `vertices_edges` property.
    /// However, note that the `vertices_edges` arrays will not be sorted in counterclockwise order.
    pub fn edges_vertices_to_vertices_edges_unsorted(&mut self) {
        let num_vertices = self.num_vertices();
        self.vertices_edges = Some(edges_vertices_to_vertices_edges_unsorted(
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
            num_vertices
        ));
    }

    /// Given a FOLD object with `vertices_coords_f64` or `vertices_coords_exact` is `Some(2D vectors)`
    /// and `vertices_vertices == Some(_)`,
    /// sorts each vertices_vertices array in counterclockwise order around the vertex in the plane,
    /// according to `AngleRep::angle_rep`.
    pub fn sort_vertices_vertices(&mut self) {
        self.vertices_vertices = Some(sort_vertices_vertices(
            self.vertices_vertices.take().expect("vertices_vertices must exist"),
            self.coords_ref().expect("vertices_coords_exact or vertices_coords_f64 must exist")
        ));
    }

    /// Given a FOLD object frame with `vertices_coords_f64` or `vertices_coords_exact` is `Some(2D vectors)`
    /// and `edges_vertices == Some(_)` (defining edge endpoints),
    /// automatically computes the `vertices_vertices` property and
    /// sorts them counterclockwise by angle in the plane.
    /// according to `AngleRep::angle_rep()`.
    pub fn edges_vertices_to_vertices_vertices_sorted(&mut self) {
        self.vertices_vertices = Some(edges_vertices_to_vertices_vertices_sorted(
            self.coords_ref().expect("vertices_coords_exact or vertices_coords_f64 must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
        ));
    }

    /// Given a FOLD object frame with `vertices_coords_f64` or `vertices_coords_exact` is `Some(2D vectors)`
    /// and `edges_vertices == Some(_)` (defining edge endpoints),
    /// automatically computes the `vertices_vertices` and `vertices_edges` properties and
    /// sorts them counterclockwise by angle in the plane.
    /// according to `AngleRep::angle_rep()`.
    pub fn edges_vertices_to_vertices_edges_sorted(&mut self) {
        let (vertices_vertices, vertices_edges) = edges_vertices_to_vertices_edges_sorted(
            self.coords_ref().expect("vertices_coords_exact or vertices_coords_f64 must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
        );
        self.vertices_vertices = Some(vertices_vertices);
        self.vertices_edges = Some(vertices_edges);
    }

    /// Given a FOLD object frame with `vertices_vertices == Some(_)` and `edges_vertices == Some(_)`,
    /// fills in the corresponding `vertices_edges` property (preserving order).
    pub fn vertices_vertices_to_vertices_edges(&mut self) {
        self.vertices_edges = Some(vertices_vertices_to_vertices_edges(
            self.vertices_vertices.as_ref().expect("vertices_vertices must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
        ));
    }

    /// Given a FOLD object frame with counterclockwise-sorted `vertices_vertices`,
    /// constructs the implicitly defined faces, setting the `faces_vertices` property.
    /// If `vertices_coords_f64` or `vertices_coords_exact` is `Some(2D points)`, then
    /// this excludes faces that are not counterclockwise.
    pub fn vertices_vertices_to_faces_vertices(&mut self) {
        self.faces_vertices = Some(vertices_vertices_to_faces_vertices(
            self.vertices_vertices.as_ref().expect("vertices_vertices must exist"),
            self.coords_ref()
        ));
    }

    /// Given a FOLD object frame with counterclockwise-sorted `vertices_edges` property,
    /// **and with `edges_vertices = Some(_)`**,
    /// constructs the implicitly defined faces, setting both `faces_vertices` and `faces_edges` properties.
    /// Handles multiple edges between the same two vertices (unlike `convert::vertices_vertices_to_faces_vertices`).
    /// If `vertices_coords_f64` or `vertices_coords_exact` is `Some(2D points)`, then
    /// this excludes faces that are not counterclockwise.
    pub fn vertices_edges_to_faces_vertices_edges(&mut self) {
        let (faces_vertices, faces_edges) = vertices_edges_to_faces_vertices_edges(
            self.vertices_edges.as_ref().expect("vertices_edges must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
            self.coords_ref()
        );
        self.faces_vertices = Some(faces_vertices);
        self.faces_edges = Some(faces_edges);
    }

    /// Given a FOLD object frame with `vertices_coords_exact` or `vertices_coords_f64` is `Some(2D points)`
    /// and `edges_vertices == Some(_)`,
    /// computes a counterclockwise-sorted `vertices_vertices` property
    /// and constructs the implicitly defined faces, setting `faces_vertices` property.
    pub fn edges_vertices_to_faces_vertices(&mut self) {
        let (vertices_vertices, faces_vertices) = edges_vertices_to_faces_vertices(
            self.coords_ref().expect("vertices_coords_exact or vertices_coords_f64 must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
        );
        self.vertices_vertices = Some(vertices_vertices);
        self.faces_vertices = Some(faces_vertices);
    }

    /// Given a FOLD object frame with `vertices_coords_exact` or `vertices_coords_f64` is `Some(2D points)`
    /// and `edges_vertices == Some(_)`,
    /// computes counterclockwise-sorted `vertices_vertices` and `vertices_edges` properties and
    /// constructs the implicitly defined faces, setting both `faces_vertices` and `faces_edges` property.
    pub fn edges_vertices_to_faces_vertices_edges(&mut self) {
        let (vertices_vertices, vertices_edges, faces_vertices, faces_edges) = edges_vertices_to_faces_vertices_edges(
            self.coords_ref().expect("vertices_coords_exact or vertices_coords_f64 must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
        );
        self.vertices_vertices = Some(vertices_vertices);
        self.vertices_edges = Some(vertices_edges);
        self.faces_vertices = Some(faces_vertices);
        self.faces_edges = Some(faces_edges);
    }
    
    /// Given a FOLD object frame with `faces_vertices == Some(_)`
    /// and `edges_vertices == Some(_)`,
    /// fills in the corresponding faces_edges property (preserving order).
    pub fn faces_vertices_to_faces_edges(&mut self) {
        self.faces_edges = Some(faces_vertices_to_faces_edges(
            self.faces_vertices.as_ref().expect("faces_vertices must exist"),
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
        ));
    }

    /// Given a FOLD object frame with `faces_vertices == Some(_)`, automatically computes
    /// `edges_vertices`, `edges_faces`, `faces_edges`, and `edges_assignment`.
    /// (indicating which edges are boundary with `EdgeAssignment::Boundary` and marking
    /// all other edges with `EdgeAssignment::Unassigned`).
    /// Currently assumes an orientable manifold, and uses `None`s to
    /// represent missing neighbor faces in `edges_faces` (for boundary edges).
    pub fn faces_vertices_to_edges(&mut self) {
        let (edges_vertices, edges_faces, faces_edges, edges_assignment) = faces_vertices_to_edges(
            self.faces_vertices.as_ref().expect("faces_vertices must exist"),
        );
        self.edges_vertices = Some(edges_vertices);
        self.edges_faces = Some(edges_faces);
        self.faces_edges = Some(faces_edges);
        self.edges_assignment = Some(edges_assignment);
    }

    /// Given a FOLD object frame with `edges_vertices == Some(_)`
    /// and `faces_vertices == Some(_)`,
    /// fills in `edges_faces` and `faces_edges`.
    /// Currently assumes an orientable manifold, and uses `None`s to
    /// represent missing neighbor faces in `edges_faces` (for boundary edges).
    pub fn edges_vertices_to_edges_faces_edges(&mut self) {
        let (edges_faces, faces_edges) = edges_vertices_to_edges_faces_edges(
            self.edges_vertices.as_ref().expect("edges_vertices must exist"),
            self.faces_vertices.as_ref().expect("faces_vertices must exist"),
        );
        self.edges_faces = Some(edges_faces);
        self.faces_edges = Some(faces_edges);
    }
}

#[cfg(test)]
mod test {
    use exact_number::based_expr;
    use nalgebra::DMatrix;

    use crate::fold::Frame;

    #[test]
    fn test_edges_vertices_to_vertices_vertices_no_edges() {
        let mut frame = Frame {
            edges_vertices: Some(vec![]),
            ..Default::default()
        };
        frame.edges_vertices_to_vertices_vertices_unsorted();
        assert_eq!(frame.vertices_vertices, Some(vec![]));
    }

    #[test]
    fn test_edges_vertices_to_vertices_vertices_star() {
        let mut frame = Frame {
            edges_vertices: Some(vec![
                [0, 1],
                [0, 2],
                [0, 3],
            ]),
            ..Default::default()
        };
        frame.edges_vertices_to_vertices_vertices_unsorted();
        frame.vertices_vertices.as_mut().unwrap()[0].sort();
        assert_eq!(frame.vertices_vertices, Some(vec![
            vec![1, 2, 3],
            vec![0],
            vec![0],
            vec![0],
        ]));
    }

    #[test]
    fn test_edges_vertices_to_vertices_vertices_square_slash() {
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
        frame.edges_vertices_to_vertices_vertices_unsorted();
        frame.vertices_vertices.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(frame.vertices_vertices, Some(vec![
            vec![1, 3],
            vec![0, 2, 3],
            vec![1, 3],
            vec![0, 1, 2],
        ]));
    }

    #[test]
    fn test_edges_vertices_to_vertices_vertices_isolated_vertex() {
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
        frame.edges_vertices_to_vertices_vertices_unsorted();
        frame.vertices_vertices.as_mut().unwrap()[0].sort();
        assert_eq!(frame.vertices_vertices, Some(vec![
            vec![1, 2],
            vec![0],
            vec![0],
            vec![],
        ]));
    }

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