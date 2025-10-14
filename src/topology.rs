use std::fmt::Display;

use indexmap::IndexMap;
use serde_json::Value;

use crate::{fold::Frame, ser_de::{SerDeFold, SerDeFrame}};

/// Given `vertices_edges` (defining edge endpoints),
/// tries to compute the `edges_vertices` property.
/// The sort order of vertices in an edge is arbitrary.
pub fn try_ve_to_ev(vertices_edges: &[Vec<usize>]) -> Result<Vec<[usize; 2]>, Vec<TopologyError>> {
    let num_edges = vertices_edges.iter().flatten().copied().max().map(|n | n + 1).unwrap_or(0);
    let mut edges_vertices = vec![vec![]; num_edges];
    let mut errors = vec![];
    for (v, edges) in vertices_edges.iter().enumerate() {
        for &e in edges.iter() {
            edges_vertices[e].push(v);
        }
    }
    for (e, vertices) in edges_vertices.iter().enumerate() {
        if vertices.len() != 2 {
            errors.push(TopologyError::EdgeDoesNotHave2Vertices { edge: e, vertices: vertices.clone() });
        } else if vertices[0] == vertices[1] {
            errors.push(TopologyError::EdgeIsSelfLoop{ edge: e, vertex: vertices[0] });
        }
    }
    if errors.len() > 0 { Err(errors)?  }
    
    Ok(edges_vertices.into_iter().map(|vs| [vs[0], vs[1]]).collect::<Vec<_>>())
}

/// Given `edges_vertices` (defining edge endpoints),
/// tries to compute the `vertices_edges` property.
/// However, note that the `vertices_edges` arrays will not be sorted in counterclockwise order.
pub fn try_ev_to_ve(edges_vertices: &[[usize; 2]], num_vertices: usize) -> Result<Vec<Vec<usize>>, Vec<TopologyError>> {
    let mut vertices_edges = vec![vec![]; num_vertices];
    let mut errors = vec![];
    for (e, &[v0, v1]) in edges_vertices.iter().enumerate() {
        if v0 == v1 {
            errors.push(TopologyError::EdgeIsSelfLoop{ edge: e, vertex: v0 });
        }
        vertices_edges[v0].push(e);
        vertices_edges[v1].push(e);
    }
    if errors.is_empty() {
        Ok(vertices_edges)
    } else {
        Err(errors)
    }
}

/// Given `faces_edges` and `edges_vertices`
/// tries to calculate `faces_vertices_edges`
pub fn try_fe_ev_to_fve(faces_edges: &[Vec<usize>], edges_vertices: &[[usize; 2]])
    -> Result<Vec<Vec<(usize, usize)>>, Vec<TopologyError>>
{
    todo!()
}

/// Given `faces_vertices` and `vertices_edges`
/// tries to calculate `faces_vertices_edges`
pub fn try_fv_ve_to_fve(faces_vertices: &[Vec<usize>], vertices_edges: &[Vec<usize>])
    -> Result<Vec<Vec<(usize, usize)>>, Vec<TopologyError>>
{
    todo!()
}

/// Given `faces_vertices` and `faces_edges`
/// tries to calculate `faces_vertices_edges`
pub fn try_fv_fe_to_fve(faces_vertices: &[Vec<usize>], faces_edges: &[Vec<usize>])
    -> Result<Vec<Vec<(usize, usize)>>, Vec<TopologyError>>
{
    todo!()
}

/// Given `faces_vertices`
/// tries to calculate `edges_vertices`, making up a numbering for the edges.
pub fn try_fv_to_ev(faces_vertices: &[Vec<usize>])
    -> Result<Vec<[usize; 2]>, Vec<TopologyError>>
{
    todo!()
}

/// Given `faces_vertices_edges`
/// tries to calculate `edges_vertices`
pub fn try_fve_to_ev(faces_vertices_edges: &[Vec<(usize, usize)>])
    -> Result<Vec<[usize; 2]>, Vec<TopologyError>>
{
    todo!()
}

/// Given `faces_vertices_edges`
/// tries to calculate `edges_faces`
pub fn try_fve_to_ef(faces_vertices_edges: &[Vec<(usize, usize)>], num_edges: usize)
    -> Result<Vec<Vec<Option<usize>>>, Vec<TopologyError>>
{
    todo!()
}

/// Given `vertices_vertices`
/// tries to calculate `edges_vertices`, making up a numbering for the edges.
pub fn try_vv_to_ev(vertices_vertices: &[Vec<usize>])
    -> Result<Vec<[usize; 2]>, Vec<TopologyError>>
{
    todo!()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Element {
    Vertex,
    Edge,
    Face
}

impl Display for Element {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Vertex => write!(f, "vertex"),
            Self::Edge => write!(f, "edge"),
            Self::Face => write!(f, "face"),
        }
    }
}

/// An error in the topology
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TopologyError {
    /// The number of elements in two fields don't match.
    /// The fields, their reported counts, and whether they're just lower bounds or exact are given.
    ElementCountMismatch { element: Element, fields: [String; 2], counts: [usize; 2], exact: [bool; 2] },
    /// An edges doesn't have exactly 2 vertices. The edge and vertices are given.
    EdgeDoesNotHave2Vertices { edge: usize, vertices: Vec<usize> },
    /// An edge is a self-loop. The edge and vertex are given.
    EdgeIsSelfLoop { edge: usize, vertex: usize },
    /// Information for a vertex exists without the necessary information to support it.
    InvalidVertexInformationExistence { field: String },
    /// Information for an edge exists without the necessary information to support it.
    InvalidEdgeInformationExistence { field: String },
    /// Information for a face exists without the necessary information to support it.
    InvalidFaceInformationExistence { field: String },
}

impl Display for TopologyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ElementCountMismatch { element, fields, counts, exact } => {
                write!(f, "{element} counts aren't compatible between \"{}\" (count {} {}) and \"{}\" (count {} {})",
                    fields[0], if exact[0] { "=" } else { "≥" }, counts[0],
                    fields[1], if exact[1] { "=" } else { "≥" }, counts[1])
            },
            Self::EdgeDoesNotHave2Vertices { edge, vertices } => {
                write!(f, "edge {} does not have 2 vertices (its vertices are {:?})", edge, vertices)
            },
            Self::EdgeIsSelfLoop { edge, vertex } => {
                write!(f, "edge {} is a self-loop (vertex {})", edge, vertex)
            },
            Self::InvalidVertexInformationExistence { field } => {
                write!(f, concat!("vertex information exists (\"{}\") without anything relating vertices to edges/faces ",
                    "(\"edges_vertices\", \"vertices_edges\", or \"faces_vertices\")"), field)
            }
            Self::InvalidEdgeInformationExistence { field } => {
                write!(f, concat!("edge information exists (\"{}\") without anything relating edges to vertices/faces ",
                    "(\"edges_vertices\", \"vertices_edges\", or \"faces_edges\")"), field)
            }
            Self::InvalidFaceInformationExistence { field } => {
                write!(f, concat!("face information exists (\"{}\") without anything relating faces to vertices/edges ",
                    "(\"faces_vertices\" or \"faces_edges\")"), field)
            }
        }
    }
}

impl SerDeFrame {
    fn num_elements(fields: impl IntoIterator<Item = (String, usize, bool)>, element: Element) -> Result<Option<(usize, String)>, Vec<TopologyError>> {
        let mut errors = vec![];
        let max = fields.into_iter().reduce( |mut acc, mut elem|
            {
                if !acc.2 && elem.2 { std::mem::swap(&mut acc, &mut elem) }
                match (acc.2, elem.2) {
                    (false, false) => if acc.1 >= elem.1 { acc } else { elem }
                    (false, true) => unreachable!(),
                    (true, _) => if if elem.2 { acc.1 == elem.1 } else { acc.1 >= elem.1 } { acc } else {
                        errors.push(TopologyError::ElementCountMismatch {
                            element,
                            fields: [acc.0.clone(), elem.0.clone()],
                            counts: [acc.1, elem.1],
                            exact: [acc.2, elem.2],
                        });
                        acc
                    }
                }
            })
            .map(|(a, b, c)| (a, Some(b), c)).unwrap_or((String::new(), None, true));

        if errors.is_empty() { Ok(max.1.map(|count| (count, max.0))) } else { Err(errors) }
    }
    
    /// Counts the number of vertices
    /// Looks in fields labelled `vertices_*` (reserved and custom),
    /// then  in fields labelled `*_vertices` (reserved and custom)
    /// 
    /// Returns None if it doesn't find a vertex field.
    /// 
    /// Returns an `Err` if the fields aren't consistent about the number of vertices.
    pub(crate) fn num_vertices(&self, vertices_custom: &IndexMap<String, Vec<Value>>) -> Result<Option<(usize, String)>, Vec<TopologyError>> {
        let counts = vec![
            self.vertices_coords_f64  .as_ref().map(|v| ("vertices_coords".to_owned(),       v.len(), true)),
            self.vertices_coords_exact.as_ref().map(|v| ("vertices_exact:coords".to_owned(), v.len(), true)),
            self.vertices_vertices    .as_ref().map(|v| ("vertices_vertices".to_owned(),     v.len(), true)),
            self.vertices_edges       .as_ref().map(|v| ("vertices_edges".to_owned(),        v.len(), true)),
            self.vertices_faces       .as_ref().map(|v| ("vertices_faces".to_owned(),        v.len(), true)),
            self.edges_vertices.as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| ("edges_vertices".to_owned(), n + 1, false)),
            self.faces_vertices.as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| ("faces_vertices".to_owned(), n + 1, false)),
        ].into_iter().flatten().chain(
            vertices_custom.iter().map(|(k, v)| (format!("vertices_{k}"), v.len(), true))
        );
        Self::num_elements(counts, Element::Vertex)
    }

    /// Counts the number of edges
    /// Looks in fields labelled `edges_*` (reserved and custom),
    /// then  in fields labelled `*_edges` (reserved and custom)
    /// 
    /// Returns None if it doesn't find an edge field.
    /// 
    /// Returns an `Err` if the fields aren't consistent about the number of edges.
    pub(crate) fn num_edges(&self, edges_custom: &IndexMap<String, Vec<Value>>) -> Result<Option<(usize, String)>, Vec<TopologyError>> {
        let counts = vec![
            self.edges_vertices      .as_ref().map(|v| ("edges_vertices"       .to_owned(), v.len(), true)),
            self.edges_faces         .as_ref().map(|v| ("edges_faces"          .to_owned(), v.len(), true)),
            self.edges_assignment    .as_ref().map(|v| ("edges_assignment"     .to_owned(), v.len(), true)),
            self.edges_fold_angle_f64.as_ref().map(|v| ("edges_foldAngle"      .to_owned(), v.len(), true)),
            self.edges_fold_angle_cs .as_ref().map(|v| ("edges_exact:foldAngle".to_owned(), v.len(), true)),
            self.edges_length_f64    .as_ref().map(|v| ("edges_length"         .to_owned(), v.len(), true)),
            self.edges_length2_exact .as_ref().map(|v| ("edges_exact:length2"  .to_owned(), v.len(), true)),
            self.vertices_edges.as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| ("vertices_edges".to_owned(), n + 1, false)),
            self.faces_edges   .as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| ("faces_edges"   .to_owned(), n + 1, false)),
            self.edge_orders   .as_ref().and_then(|v| v.iter().flat_map(|&(a, b, _)| [a, b]).max()).map(|n| ("edge_orders".to_owned(), n + 1, false)),
        ].into_iter().flatten().chain(
            edges_custom.iter().map(|(k, v)| (format!("edges_{k}"), v.len(), true))
        );
        Self::num_elements(counts, Element::Edge)
    }

    /// Counts the number of faces
    /// Looks in fields labelled `faces_*` (reserved and custom),
    /// then  in fields labelled `*_faces` (reserved and custom)
    /// 
    /// Returns None if it doesn't find an face field.
    /// 
    /// Returns an `Err` if the fields aren't consistent about the number of faces.
    pub(crate) fn num_faces(&self, faces_custom: &IndexMap<String, Vec<Value>>) -> Result<Option<(usize, String)>, Vec<TopologyError>> {
        let counts = vec![
            self.faces_vertices.as_ref().map(|v| ("faces_vertices".to_owned(), v.len(), true)),
            self.faces_edges   .as_ref().map(|v| ("faces_edges"   .to_owned(), v.len(), true)),
            self.faces_faces   .as_ref().map(|v| ("faces_faces"   .to_owned(), v.len(), true)),
            self.vertices_faces.as_ref().and_then(|v| v.iter().flatten().flatten().copied().max()).map(|n| ("vertices_faces".to_owned(), n + 1, false)),
            self.edges_faces   .as_ref().and_then(|v| v.iter().flatten().flatten().copied().max()).map(|n| ("edges_faces"   .to_owned(), n + 1, false)),
            self.face_orders.as_ref().and_then(|v| v.iter().flat_map(|&(a, b, _)| [a, b]).max()).map(|n| ("faces_orders".to_owned(), n + 1, false)),
        ].into_iter().flatten().chain(
            faces_custom.iter().map(|(k, v)| (format!("faces_{k}"), v.len(), true))
        );
        Self::num_elements(counts, Element::Face)
    }

    /// Extracts `vertices_edges`, `edges_vertices`, `edges_faces`, and `faces_vertices_edges` from this intermediate type.
    pub(crate) fn extract_topology(mut self, vertices_meta: Option<(usize, String)>, edges_meta: Option<(usize, String)>, faces_meta: Option<(usize, String)>)
        -> Result<(Self, Option<Vec<Vec<usize>>>, Option<Vec<[usize; 2]>>, Option<Vec<Vec<Option<usize>>>>, Option<Vec<Vec<(usize, usize)>>>), Vec<TopologyError>>
    {
        let num_vertices = vertices_meta.as_ref().map(|&(n, _)| n).unwrap_or(0);
        
        let (edges_vertices, vertices_edges) = match (self.edges_vertices.take(), self.vertices_edges.take()) {
            (None, None) => (None, None),
            (None, Some(ve)) => (Some(try_ve_to_ev(&ve)?), Some(ve)),
            (Some(ev), _) => { let t = (Some(try_ev_to_ve(&ev, num_vertices)?), Some(ev)); (t.1, t.0) }
        };

        let faces_vertices = self.faces_vertices.take();
        let faces_edges = self.faces_edges.take();

        let (vertices_edges, edges_vertices, edges_faces, faces_vertices_edges) =
            if let (Some(vertices_edges), Some(edges_vertices)) = (vertices_edges, edges_vertices) {
                if let Some(faces_edges) = faces_edges {
                    let fve = try_fe_ev_to_fve(&faces_edges, &edges_vertices)?;
                    let ef = try_fve_to_ef(&fve, edges_vertices.len())?;
                    (Some(vertices_edges), Some(edges_vertices), Some(ef), Some(fve))

                } else if let Some(faces_vertices) = faces_vertices {
                    let fve = try_fv_ve_to_fve(&faces_vertices, &vertices_edges)?;
                    let ef = try_fve_to_ef(&fve, edges_vertices.len())?;
                    (Some(vertices_edges), Some(edges_vertices), Some(ef), Some(fve))

                } else {
                    // Linkage;  no face information!
                    if let Some((_, field)) = faces_meta {
                        Err(vec![TopologyError::InvalidFaceInformationExistence { field }])?;
                    }
                    (Some(vertices_edges), Some(edges_vertices), None, None)
                }
            } else {
                if let Some(faces_vertices) = faces_vertices {
                    if let Some(faces_edges) = faces_edges {
                        let fve = try_fv_fe_to_fve(&faces_vertices, &faces_edges)?;
                        let ev = try_fve_to_ev(&fve)?;
                        let ve = try_ev_to_ve(&ev, num_vertices)?;
                        let ef = try_fve_to_ef(&fve, ev.len())?;
                        (Some(ve), Some(ev), Some(ef), Some(fve))
                    } else {
                        // No edge information!
                        if let Some((_, field)) = edges_meta {
                            Err(vec![TopologyError::InvalidEdgeInformationExistence { field }])?;
                        }
                        let ev = try_fv_to_ev(&faces_vertices)?;
                        let ve = try_ev_to_ve(&ev, num_vertices)?;
                        let fve = try_fv_ve_to_fve(&faces_vertices, &ve)?;
                        let ef = try_fve_to_ef(&fve, ev.len())?;
                        (Some(ve), Some(ev), Some(ef), Some(fve))
                    }
                } else if let Some(vertices_vertices) = self.vertices_vertices.take() {
                    // No edge or face information!
                    let mut errors = vec![];
                    if let Some((_, field)) = edges_meta { errors.push(TopologyError::InvalidEdgeInformationExistence { field }); }
                    if let Some((_, field)) = faces_meta { errors.push(TopologyError::InvalidFaceInformationExistence { field }); }
                    if !errors.is_empty() { Err(errors)? }
                    let ev = try_vv_to_ev(&vertices_vertices)?;
                    let ve = try_ev_to_ve(&ev, num_vertices)?;
                    (Some(ve), Some(ev), None, None)
                } else {
                    let mut errors = vec![];
                    if let Some((_, field)) = vertices_meta { errors.push(TopologyError::InvalidVertexInformationExistence { field }); }
                    if let Some((_, field)) = edges_meta { errors.push(TopologyError::InvalidEdgeInformationExistence { field }); }
                    if let Some((_, field)) = faces_meta { errors.push(TopologyError::InvalidFaceInformationExistence { field }); }
                    if !errors.is_empty() { Err(errors)? }
                    (None, None, None, None)
                }
            };
            
        Ok((self, vertices_edges, edges_vertices, edges_faces, faces_vertices_edges))
    }
}

impl Frame {
    /// Extracts `vertices_vertices`, `vertices_edges`, `vertices_faces`,
    /// `edges_vertices`, `edges_faces`,
    /// `faces_vertices`, `faces_edges`, and `faces_faces` from this frame.
    pub fn extract_topology(self)
        -> Result<(Self, Option<Vec<Vec<usize>>>, Option<Vec<Vec<usize>>>, Option<Vec<Vec<Option<usize>>>>,
                         Option<Vec<[usize; 2]>>,                          Option<Vec<Vec<Option<usize>>>>,
                         Option<Vec<Vec<usize>>>, Option<Vec<Vec<usize>>>, Option<Vec<Vec<Option<usize>>>>), String>
    {
        unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use exact_number::{based_expr, rat::Rat};
    use indexmap::IndexMap;
    use nalgebra::DMatrix;
    use serde_json::json;

    use crate::{filter::{edges_vertices_incident, other_vertex}, ser_de::{SerDeFrame, SerDeRat}, topology::TopologyError};

    #[test]
    fn test_num_vertices_vertices_coords_f64() {
        let frame = SerDeFrame {
            vertices_coords_f64: Some(vec![
                vec![1.0, 2.0],
                vec![3.0, 4.0],
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((2, "vertices_coords".to_owned()))));
    }

    #[test]
    fn test_num_vertices_vertices_coords_exact() {
        let frame = SerDeFrame {
            vertices_coords_exact: Some(vec![
                vec![vec![SerDeRat(Rat::ONE)], vec![SerDeRat(Rat::TWO)]],
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((1, "vertices_exact:coords".to_owned()))));
    }

    #[test]
    fn test_num_vertices_vertices_vertices() {
        let frame = SerDeFrame {
            vertices_vertices: Some(vec![
                vec![1, 4, 5],
                vec![3, 4, 1, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((2, "vertices_vertices".to_owned()))));
    }

    #[test]
    fn test_num_vertices_vertices_edges() {
        let frame = SerDeFrame {
            vertices_edges: Some(vec![
                vec![1, 4, 5],
                vec![3, 4, 1, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((2, "vertices_edges".to_owned()))));
    }

    #[test]
    fn test_num_vertices_vertices_edges_mismatched() {
        let frame = SerDeFrame {
            vertices_vertices: Some(vec![
                vec![1, 2],
            ]),
            vertices_edges: Some(vec![
                vec![1, 4, 5],
                vec![3, 4, 1, 2]
            ]),
            ..Default::default()
        };
        assert!(frame.num_vertices(&IndexMap::new()).is_err());
    }

    #[test]
    fn test_num_vertices_vertices_faces() {
        let frame = SerDeFrame {
            vertices_faces: Some(vec![
                vec![Some(1), None, Some(3)],
                vec![Some(0), Some(3), Some(4)],
                vec![None, Some(2), Some(5)],
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((3, "vertices_faces".to_owned()))));
    }

    #[test]
    fn test_num_vertices_vertices_custom() {
        let frame = SerDeFrame {
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&vec![
                ("test:field".to_owned(), vec![
                    json!(3),
                    json!(5),
                    json!(5),
                ])
            ].into_iter().collect::<IndexMap<_, _>>()), Ok(Some((3, "vertices_test:field".to_owned()))));
    }

    #[test]
    fn test_num_vertices_edges_vertices() {
        let frame = SerDeFrame {
            edges_vertices: Some(vec![
                [3, 4],
                [1, 4],
                [5, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((6, "edges_vertices".to_owned()))));
    }

    #[test]
    fn test_num_vertices_faces_vertices() {
        let frame = SerDeFrame {
            faces_vertices: Some(vec![
                vec![0, 3, 4],
                vec![2, 1, 4],
                vec![3, 5, 2]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((6, "faces_vertices".to_owned()))));
    }

    #[test]
    fn test_num_vertices_edges_vertices_faces_vertices() {
        let mut frame = SerDeFrame {
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
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((8, "faces_vertices".to_owned()))));
        frame.edges_vertices = Some(vec![
            [10, 2],
            [3, 5]
        ]);
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((11, "edges_vertices".to_owned()))));
    }

    #[test]
    fn test_num_vertices_no_vertex_field() {
        let frame = SerDeFrame {
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(None));
    }

    #[test]
    fn test_try_ev_to_ve() {
        let vertices_edges = super::try_ev_to_ve(&[], 0);
        assert_eq!(vertices_edges, Ok(vec![]));

        // Star
        let mut vertices_edges = super::try_ev_to_ve(&[
            [0, 1],
            [0, 2],
            [0, 3],
        ], 4);
        vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(vertices_edges, Ok(vec![
            vec![0, 1, 2],
            vec![0],
            vec![1],
            vec![2]
        ]));

        // Slash
        let mut vertices_edges = super::try_ev_to_ve(&[
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [1, 3],
        ], 4);
        vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(vertices_edges, Ok(vec![
            vec![0, 3],
            vec![0, 1, 4],
            vec![1, 2],
            vec![2, 3, 4],
        ]));

        // Doubled edge
        let mut vertices_edges = super::try_ev_to_ve(&[
            [0, 1],
            [1, 2],
            [1, 2],
        ], 3);
        vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(vertices_edges, Ok(vec![
            vec![0],
            vec![0, 1, 2],
            vec![1, 2],
        ]));

        // Isolated vertex
        let mut vertices_edges = super::try_ev_to_ve(&[
            [0, 1],
            [0, 2],
        ], 4);
        vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(vertices_edges, Ok(vec![
            vec![0, 1],
            vec![0],
            vec![1],
            vec![],
        ]));

        // Self-loop
        let vertices_edges = super::try_ev_to_ve(&[
            [0, 1],
            [1, 2],
            [2, 2],
        ], 3);
        assert!(vertices_edges.is_err());
    }

    #[test]
    fn test_try_ve_to_ev() {
        use TopologyError::*;
        let edges_vertices = super::try_ve_to_ev(&[]);
        assert_eq!(edges_vertices, Ok(vec![]));

        // Star
        let mut edges_vertices = super::try_ve_to_ev(&[
            vec![0, 1, 2],
            vec![0],
            vec![1],
            vec![2]
        ]);
        edges_vertices.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(edges_vertices, Ok(vec![
            [0, 1],
            [0, 2],
            [0, 3],
        ]));

        // Slash
        let mut edges_vertices = super::try_ve_to_ev(&[
            vec![0, 3],
            vec![0, 1, 4],
            vec![1, 2],
            vec![2, 3, 4],
        ]);
        edges_vertices.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(edges_vertices, Ok(vec![
            [0, 1],
            [1, 2],
            [2, 3],
            [0, 3],
            [1, 3],
        ]));

        // Doubled edge
        let mut edges_vertices = super::try_ve_to_ev(&[
            vec![0],
            vec![0, 1, 2],
            vec![1, 2],
        ]);
        edges_vertices.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(edges_vertices, Ok(vec![
            [0, 1],
            [1, 2],
            [1, 2],
        ]));

        // Isolated vertex
        let mut edges_vertices = super::try_ve_to_ev(&[
            vec![0, 1],
            vec![0],
            vec![1],
            vec![],
        ]);
        edges_vertices.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(edges_vertices, Ok(vec![
            [0, 1],
            [0, 2],
        ]));

        // Self-loop
        let edges_vertices = super::try_ve_to_ev(&[
            vec![0],
            vec![0, 2],
            vec![2, 1, 1],
        ]);
        assert_eq!(edges_vertices, Err(vec![EdgeIsSelfLoop { edge: 1, vertex: 2 }]));

        // Edge index skipped
        let edges_vertices = super::try_ve_to_ev(&[
            vec![0, 1],
            vec![1, 3],
            vec![3, 0],
        ]);
        assert_eq!(edges_vertices, Err(vec![EdgeDoesNotHave2Vertices { edge: 2, vertices: vec![] }]));

        // Edge missing vertices
        let edges_vertices = super::try_ve_to_ev(&[
            vec![0, 1],
            vec![2],
            vec![2, 0],
        ]);
        assert_eq!(edges_vertices, Err(vec![EdgeDoesNotHave2Vertices { edge: 1, vertices: vec![0] }]));

        // Edge has too many vertices
        let mut edges_vertices = super::try_ve_to_ev(&[
            vec![0, 1],
            vec![1, 2],
            vec![2, 0],
            vec![1],
        ]);
        let _ =
            edges_vertices.as_mut().map_err(|e| if let EdgeDoesNotHave2Vertices{ vertices, ..} = &mut e[0] { vertices.sort() } );
        assert_eq!(edges_vertices, Err(vec![EdgeDoesNotHave2Vertices { edge: 1, vertices: vec![0, 1, 3] }]));
    }
}