use std::fmt::Display;

use indexmap::IndexMap;
use serde_json::Value;
use typed_index_collections::{ti_vec, TiSlice, TiVec};

use crate::{filter, fold::{Edge, Face, Frame, HalfEdge, Vertex}, ser_de::{SerDeFold, SerDeFrame}};

/// Given `vertices_edges` (defining edge endpoints),
/// tries to compute the `vertices_half_edges` and `edges_vertices` properties.
pub fn try_ve_to_vh_ev(vertices_edges: &TiSlice<Vertex, Vec<Edge>>) -> Result<(TiVec<Vertex, Vec<HalfEdge>>, TiVec<Edge, [Vertex; 2]>), Vec<TopologyError>> {
    let num_edges = vertices_edges.iter().flatten().copied().max().map(|n | n.0 + 1).unwrap_or(0);
    let mut vertices_half_edges = ti_vec![vec![]; vertices_edges.len()];
    let mut edges_vertices = ti_vec![vec![]; num_edges];
    let mut errors = vec![];
    for (v, edges) in vertices_edges.iter_enumerated() {
        for &e in edges.iter() {
            vertices_half_edges[v].push(HalfEdge::new(e, !edges_vertices[e].is_empty()));
            edges_vertices[e].push(v);
        }
    }
    for (e, vertices) in edges_vertices.iter_enumerated() {
        if vertices.len() != 2 {
            errors.push(TopologyError::EdgeDoesNotHave2Vertices { edge: e, vertices: vertices.clone() });
        } else if vertices[0] == vertices[1] {
            errors.push(TopologyError::EdgeIsSelfLoop{ edge: e, vertex: vertices[0] });
        }
    }

    if errors.is_empty() {
        Ok((vertices_half_edges, edges_vertices.into_iter().map(|vs| [vs[0], vs[1]]).collect::<TiVec<_, _>>()))
    } else { Err(errors) }
}

/// Given `edges_vertices` (defining edge endpoints),
/// tries to compute the `vertices_half_edges` property.
/// However, note that the `vertices_half_edges` arrays will not be sorted in counterclockwise order.
pub fn try_ev_to_vh(edges_vertices: &TiSlice<Edge, [Vertex; 2]>, num_vertices: usize) -> Result<TiVec<Vertex, Vec<HalfEdge>>, Vec<TopologyError>> {
    let mut vertices_half_edges = ti_vec![vec![]; num_vertices];
    let mut errors = vec![];
    for (e, &[v0, v1]) in edges_vertices.iter_enumerated() {
        if v0 == v1 {
            errors.push(TopologyError::EdgeIsSelfLoop{ edge: e, vertex: v0 });
        }
        vertices_half_edges[v0].push(HalfEdge::new(e, false));
        vertices_half_edges[v1].push(HalfEdge::new(e, true));
    }
    if errors.is_empty() { Ok(vertices_half_edges) } else { Err(errors) }
}

/// Given `faces_edges` and `edges_vertices`
/// tries to calculate `faces_half_edges`, which lines up with `faces_edges`.
/// 
/// For now, assumes that `edges_vertices` is valid. Call `try_ev_to_vh` to check for validity.
pub fn try_fe_ev_to_fh(faces_edges: &TiSlice<Face, Vec<Edge>>, edges_vertices: &TiSlice<Edge, [Vertex; 2]>)
    -> Result<TiVec<Face, Vec<HalfEdge>>, Vec<TopologyError>>
{
    // TODO: Check ev again?
    let mut faces_half_edges = ti_vec![vec![]; faces_edges.len()];
    let mut errors = vec![];
    for (f, edges) in faces_edges.iter_enumerated() {
        if edges.is_empty() {
            errors.push(TopologyError::FaceDoesNotHaveAtLeast2Edges { face: f, edges: edges.clone() });
            continue;
        }
        let mut vertices = edges.iter().map(|&e| edges_vertices[e]).collect::<Vec<_>>();
        let mut flipped = vec![false; edges.len()];
        // Find pair of edges that share exactly 1 vertex to start
        let start = vertices.iter().zip(vertices.iter().cycle().skip(1)).enumerate()
            .find(|(_, (vs0, vs1))| {
                // Need to preserve vertex order in original list
                let mut vs0 = **vs0; vs0.sort();
                let mut vs1 = **vs1; vs1.sort();
                vs0 != vs1
            })
            .map(|(e, _)| e)
            .unwrap_or(0);

        // Align and chain the edges
        if filter::edges_vertices_incident(vertices[start], vertices[(start + 1) % edges.len()]) == Some(vertices[start][0]) {
            flipped[start] = !flipped[start];
            vertices[start].reverse();
        }
        for (e0, e1) in (0..edges.len()).map(|e| ((start + e) % edges.len(), (start + e + 1) % edges.len())) {
            if let Some(other) = filter::try_other_vertex(vertices[e1], vertices[e0][1]) {
                if vertices[e1][0] == other {
                    flipped[e1] = !flipped[e1];
                    vertices[e1].reverse()
                }
            } else {
                errors.push(TopologyError::ConsectiveEdgesOfFaceDoNotShareVertex {
                    face: f, edges: [edges[e0], edges[e1]], vertices: [vertices[e0], vertices[e1]]
                });
            }
        }

        faces_half_edges[f] = edges.iter().enumerate().map(|(i, &e)| HalfEdge::new(e, flipped[i])).collect::<Vec<_>>();
    }
    if errors.is_empty() { Ok(faces_half_edges) } else { Err(errors) }
}

/// Given `faces_vertices` and `vertices_half_edges`
/// tries to calculate `faces_half_edges`
/// 
/// For now, assumes that `vertices_half_edges` is valid.
pub fn try_fv_vh_to_fh(faces_vertices: &TiSlice<Face, Vec<Vertex>>, vertices_half_edges: &TiSlice<Vertex, Vec<HalfEdge>>)
    -> Result<TiVec<Face, Vec<HalfEdge>>, Vec<TopologyError>>
{
    let mut faces_half_edges = ti_vec![vec![]; faces_vertices.len()];
    let mut errors = vec![];
    for (f, vertices) in faces_vertices.iter_enumerated() {
        if vertices.is_empty() {
            errors.push(TopologyError::FaceDoesNotHaveAtLeast2Vertices { face: f, vertices: vertices.clone() });
            continue;
        }
        for (&v0, &v1) in vertices.iter().zip(vertices.iter().cycle().skip(1)) {
            let v_half_edges = filter::vertices_half_edges_incident(&vertices_half_edges[v0], &vertices_half_edges[v1]);
            if v_half_edges.len() == 0 {
                errors.push(TopologyError::ConsectiveVerticesOfFaceDoNotShareEdge { face: f, vertices: [v0, v1] })
            } else if v_half_edges.len() == 1 {
                faces_half_edges[f].push(v_half_edges[0]);
            } else {
                errors.push(TopologyError::AmbiguousDoubleEdge { face: f, edges: v_half_edges.into_iter().map(|h| h.edge()).collect(), vertices: [v0, v1] })
            }
        }
    }

    if errors.is_empty() { Ok(faces_half_edges) } else { Err(errors) }
}

/// Given `faces_vertices` and `faces_edges`
/// tries to calculate `edges_vertices` and `faces_half_edges`
/// 
/// Assumes the face counts actually match.
pub fn try_fv_fe_to_ev_fh(faces_vertices: &TiSlice<Face, Vec<Vertex>>, faces_edges: &TiSlice<Face, Vec<Edge>>, num_edges: usize)
    -> Result<(TiVec<Edge, [Vertex; 2]>, TiVec<Face, Vec<HalfEdge>>), Vec<TopologyError>>
{
    let mut errors = vec![];
    let mut edges_vertices = ti_vec![[Vertex(usize::MAX); 2]; num_edges];
    let mut faces_half_edges = ti_vec![vec![]; faces_edges.len()];

    for ((f, vertices), edges) in faces_vertices.iter_enumerated().zip(faces_edges.iter()) {
        if vertices.len() != edges.len() {
            errors.push(TopologyError::FaceVertexEdgeCountMismatch { face: f, vertices: vertices.clone(), edges: edges.clone() });
        } else if vertices.len() < 2 {
            errors.push(TopologyError::FaceDoesNotHaveAtLeast2Vertices { face: f, vertices: vertices.clone() })
        }
        for ((&v0, &e), &v1) in vertices.iter().zip(edges.iter()).zip(vertices.iter().cycle().skip(1)) {
            if edges_vertices[e][0] == Vertex(usize::MAX) {
                edges_vertices[e] = [v0, v1];
                faces_half_edges[f].push(HalfEdge::new(e, false));
            } else if edges_vertices[e] == [v0, v1] {
                faces_half_edges[f].push(HalfEdge::new(e, false));
            } else if edges_vertices[e] == [v1, v0] {
                faces_half_edges[f].push(HalfEdge::new(e, true));
            } else {
                errors.push(TopologyError::EdgeVerticesMismatch { face: f, edge: e, vertices: [edges_vertices[e], [v0, v1]] })
            }
        }
    }

    let mut unaccounted_edges = edges_vertices.iter_enumerated()
        .filter(|&(_, vs)| vs[0] == Vertex(usize::MAX))
        .map(|(e, _)| e)
        .peekable();
    if unaccounted_edges.peek().is_some() {
        errors.push(TopologyError::AmbiguousMissingEdgesFromFaces { missing_edges: unaccounted_edges.collect() })
    }

    if errors.is_empty() { Ok((edges_vertices, faces_half_edges)) } else { Err(errors) }
}

/// Given `faces_vertices`
/// tries to calculate `edges_vertices`, making up a numbering for the edges.
pub fn try_fv_to_ev(faces_vertices: &TiSlice<Face, Vec<Vertex>>)
    -> Result<TiVec<Edge, [Vertex; 2]>, Vec<TopologyError>>
{
    todo!()
}

/// Given `faces_half_edges`
/// tries to calculate `edges_vertices`
pub fn try_fh_to_ev(faces_half_edges: &TiSlice<Face, Vec<HalfEdge>>)
    -> Result<TiVec<Edge, [Vertex; 2]>, Vec<TopologyError>>
{
    todo!()
}

/// Given `faces_half_edges`
/// tries to calculate `edges_faces`
pub fn try_fh_to_ef(faces_half_edges: &TiSlice<Face, Vec<HalfEdge>>, num_edges: usize)
    -> Result<TiVec<Edge, Vec<Option<Face>>>, Vec<TopologyError>>
{
    todo!()
}

/// Given `vertices_vertices`
/// tries to calculate `edges_vertices`, making up a numbering for the edges.
pub fn try_vv_to_ev(vertices_vertices: &TiSlice<Vertex, Vec<Vertex>>)
    -> Result<TiVec<Edge, [Vertex; 2]>, Vec<TopologyError>>
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
    /// An edge doesn't have exactly 2 vertices. The edge and vertices are given.
    EdgeDoesNotHave2Vertices { edge: Edge, vertices: Vec<Vertex> },
    /// An edge is a self-loop. The edge and vertex are given.
    EdgeIsSelfLoop { edge: Edge, vertex: Vertex },
    /// A face doesn't have at least 2 edges. The face and edges are given.
    FaceDoesNotHaveAtLeast2Edges { face: Face, edges: Vec<Edge> },
    /// A face doesn't have at least 2 vertices. The face and vertices are given.
    FaceDoesNotHaveAtLeast2Vertices { face: Face, vertices: Vec<Vertex> },
    /// A face's number of vertices and number of edges don't match.
    FaceVertexEdgeCountMismatch { face: Face, vertices: Vec<Vertex>, edges: Vec<Edge> },
    /// Consecutive edges of a face do not share a vertex.
    ConsectiveEdgesOfFaceDoNotShareVertex { face: Face, edges: [Edge; 2], vertices: [[Vertex; 2]; 2] },
    /// Consecutive vertices of a face do not share an edge.
    ConsectiveVerticesOfFaceDoNotShareEdge { face: Face, vertices: [Vertex; 2] },
    /// Encountered mismatching definitions of an edge when iterating over faces.
    EdgeVerticesMismatch { face: Face, edge: Edge, vertices: [[Vertex; 2]; 2] },
    /// A double-edge was found without enough information to resolve which edge a face takes.
    /// The face, vertices, and edges are given.
    AmbiguousDoubleEdge { face: Face, edges: Vec<Edge>, vertices: [Vertex; 2] },
    /// `faces_edges` doesn't account for all the edges, and there isn't enough information
    /// to derive `edges_vertices`.
    AmbiguousMissingEdgesFromFaces { missing_edges: Vec<Edge> },
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
            Self::FaceDoesNotHaveAtLeast2Edges { face, edges } => {
                write!(f, "face {} does not have at least 2 edges (its edges are {:?})", face, edges)
            },
            Self::FaceDoesNotHaveAtLeast2Vertices { face, vertices } => {
                write!(f, "face {} does not have at least 2 vertices (its vertices are {:?})", face, vertices)
            },
            Self::FaceVertexEdgeCountMismatch { face, vertices, edges } => {
                write!(f, "face {face} has {} vertices ({vertices:?}) but {} edges ({edges:?}), so the counts don't match",
                    vertices.len(), edges.len())
            }
            Self::ConsectiveEdgesOfFaceDoNotShareVertex { face, edges, vertices } => {
                write!(f, "consecutive edges {edges:?} of face {face} do not share a vertex (their vertices are {:?}, {:?})",
                    vertices[0], vertices[1])
            }
            Self::ConsectiveVerticesOfFaceDoNotShareEdge { face, vertices } => {
                write!(f, "consecutive vertices {vertices:?} of face {face} do not share an edge")
            }
            Self::EdgeVerticesMismatch { face, edge, vertices } => {
                write!(f, "edge {edge} is defined twice ({:?}, {:?}), the second definition being encountered with face {face}",
                    vertices[0], vertices[1])
            }
            Self::AmbiguousDoubleEdge { face, edges, vertices } => {
                write!(f, concat!("consecutive vertices {:?} of face {} share multiple edges {:?}",
                    "without enough information (e.g. faces_edges) to resolve which edge belongs to the face"),
                    vertices, face, edges)
            }
            Self::AmbiguousMissingEdgesFromFaces { missing_edges } => {
                write!(f, concat!("\"faces_edges\" is missing edges {:?} without the necessary information
                    (e.g. \"edges_vertices\" or \"vertices_edges\") to derive \"edges_vertices\""), missing_edges)
            }
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
    pub(crate) fn num_vertices(&self, vertices_custom: &IndexMap<String, TiVec<Vertex, Value>>) -> Result<Option<(usize, String)>, Vec<TopologyError>> {
        let counts = vec![
            self.vertices_coords_f64  .as_ref().map(|v| ("vertices_coords".to_owned(),       v.len(), true)),
            self.vertices_coords_exact.as_ref().map(|v| ("vertices_exact:coords".to_owned(), v.len(), true)),
            self.vertices_vertices    .as_ref().map(|v| ("vertices_vertices".to_owned(),     v.len(), true)),
            self.vertices_edges       .as_ref().map(|v| ("vertices_edges".to_owned(),        v.len(), true)),
            self.vertices_faces       .as_ref().map(|v| ("vertices_faces".to_owned(),        v.len(), true)),
            self.edges_vertices.as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| ("edges_vertices".to_owned(), n.0 + 1, false)),
            self.faces_vertices.as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| ("faces_vertices".to_owned(), n.0 + 1, false)),
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
    pub(crate) fn num_edges(&self, edges_custom: &IndexMap<String, TiVec<Edge, Value>>) -> Result<Option<(usize, String)>, Vec<TopologyError>> {
        let counts = vec![
            self.edges_vertices      .as_ref().map(|v| ("edges_vertices"       .to_owned(), v.len(), true)),
            self.edges_faces         .as_ref().map(|v| ("edges_faces"          .to_owned(), v.len(), true)),
            self.edges_assignment    .as_ref().map(|v| ("edges_assignment"     .to_owned(), v.len(), true)),
            self.edges_fold_angle_f64.as_ref().map(|v| ("edges_foldAngle"      .to_owned(), v.len(), true)),
            self.edges_fold_angle_cs .as_ref().map(|v| ("edges_exact:foldAngle".to_owned(), v.len(), true)),
            self.edges_length_f64    .as_ref().map(|v| ("edges_length"         .to_owned(), v.len(), true)),
            self.edges_length2_exact .as_ref().map(|v| ("edges_exact:length2"  .to_owned(), v.len(), true)),
            self.vertices_edges.as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| ("vertices_edges".to_owned(), n.0 + 1, false)),
            self.faces_edges   .as_ref().and_then(|v| v.iter().flatten().copied().max()).map(|n| ("faces_edges"   .to_owned(), n.0 + 1, false)),
            self.edge_orders   .as_ref().and_then(|v| v.iter().flat_map(|&(a, b, _)| [a, b]).max()).map(|n| ("edge_orders".to_owned(), n.0 + 1, false)),
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
    pub(crate) fn num_faces(&self, faces_custom: &IndexMap<String, TiVec<Face, Value>>) -> Result<Option<(usize, String)>, Vec<TopologyError>> {
        let counts = vec![
            self.faces_vertices.as_ref().map(|v| ("faces_vertices".to_owned(), v.len(), true)),
            self.faces_edges   .as_ref().map(|v| ("faces_edges"   .to_owned(), v.len(), true)),
            self.faces_faces   .as_ref().map(|v| ("faces_faces"   .to_owned(), v.len(), true)),
            self.vertices_faces.as_ref().and_then(|v| v.iter().flatten().flatten().copied().max()).map(|n| ("vertices_faces".to_owned(), n.0 + 1, false)),
            self.edges_faces   .as_ref().and_then(|v| v.iter().flatten().flatten().copied().max()).map(|n| ("edges_faces"   .to_owned(), n.0 + 1, false)),
            self.face_orders.as_ref().and_then(|v| v.iter().flat_map(|&(a, b, _)| [a, b]).max()).map(|n| ("faces_orders".to_owned(), n.0 + 1, false)),
        ].into_iter().flatten().chain(
            faces_custom.iter().map(|(k, v)| (format!("faces_{k}"), v.len(), true))
        );
        Self::num_elements(counts, Element::Face)
    }

    /// Extracts `vertices_half_edges`, `edges_vertices`, `edges_faces`, and `faces_half_edges` from this intermediate type.
    pub(crate) fn extract_topology(mut self, vertices_meta: Option<(usize, String)>, edges_meta: Option<(usize, String)>, faces_meta: Option<(usize, String)>)
        -> Result<(
            Self,
            Option<TiVec<Vertex, Vec<HalfEdge>>>,
            Option<TiVec<Edge, [Vertex; 2]>>,
            Option<TiVec<Edge, Vec<Option<Face>>>>,
            Option<TiVec<Face, Vec<HalfEdge>>>
        ), Vec<TopologyError>>
    {
        let num_vertices = vertices_meta.as_ref().map(|&(n, _)| n).unwrap_or(0);
        let num_edges = edges_meta.as_ref().map(|&(n, _)| n).unwrap_or(0);
        
        let (vertices_half_edges, edges_vertices) = match (self.edges_vertices.take(), self.vertices_edges.take()) {
            (None, None) => (None, None),
            (None, Some(ve)) => {
                let (vhe, ev) = try_ve_to_vh_ev(&ve)?;
                (Some(vhe), Some(ev))
            },
            (Some(ev), _) => (Some(try_ev_to_vh(&ev, num_vertices)?), Some(ev))
        };

        let faces_vertices = self.faces_vertices.take();
        let faces_edges = self.faces_edges.take();

        let (vertices_half_edges, edges_vertices, edges_faces, faces_half_edges) =
            if let (Some(vertices_half_edges), Some(edges_vertices)) = (vertices_half_edges, edges_vertices) {
                if let Some(faces_edges) = faces_edges {
                    let fh = try_fe_ev_to_fh(&faces_edges, &edges_vertices)?;
                    let ef = try_fh_to_ef(&fh, edges_vertices.len())?;
                    (Some(vertices_half_edges), Some(edges_vertices), Some(ef), Some(fh))

                } else if let Some(faces_vertices) = faces_vertices {
                    let fh = try_fv_vh_to_fh(&faces_vertices, &vertices_half_edges)?;
                    let ef = try_fh_to_ef(&fh, edges_vertices.len())?;
                    (Some(vertices_half_edges), Some(edges_vertices), Some(ef), Some(fh))

                } else {
                    // Linkage;  no face information!
                    if let Some((_, field)) = faces_meta {
                        Err(vec![TopologyError::InvalidFaceInformationExistence { field }])?;
                    }
                    (Some(vertices_half_edges), Some(edges_vertices), None, None)
                }
            } else {
                if let Some(faces_vertices) = faces_vertices {
                    if let Some(faces_edges) = faces_edges {
                        let (ev, fh) = try_fv_fe_to_ev_fh(&faces_vertices, &faces_edges, num_edges)?;
                        let ve = try_ev_to_vh(&ev, num_vertices)?;
                        let ef = try_fh_to_ef(&fh, ev.len())?;
                        (Some(ve), Some(ev), Some(ef), Some(fh))
                    } else {
                        // No edge information!
                        if let Some((_, field)) = edges_meta {
                            Err(vec![TopologyError::InvalidEdgeInformationExistence { field }])?;
                        }
                        let ev = try_fv_to_ev(&faces_vertices)?;
                        let ve = try_ev_to_vh(&ev, num_vertices)?;
                        let fh = try_fv_vh_to_fh(&faces_vertices, &ve)?;
                        let ef = try_fh_to_ef(&fh, ev.len())?;
                        (Some(ve), Some(ev), Some(ef), Some(fh))
                    }
                } else if let Some(vertices_vertices) = self.vertices_vertices.take() {
                    // No edge or face information!
                    let mut errors = vec![];
                    if let Some((_, field)) = edges_meta { errors.push(TopologyError::InvalidEdgeInformationExistence { field }); }
                    if let Some((_, field)) = faces_meta { errors.push(TopologyError::InvalidFaceInformationExistence { field }); }
                    if !errors.is_empty() { Err(errors)? }
                    let ev = try_vv_to_ev(&vertices_vertices)?;
                    let ve = try_ev_to_vh(&ev, num_vertices)?;
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
            
        Ok((self, vertices_half_edges, edges_vertices, edges_faces, faces_half_edges))
    }
}

impl Frame {
    /// Extracts `vertices_vertices`, `vertices_edges`, `vertices_faces`,
    /// `edges_vertices`, `edges_faces`,
    /// `faces_vertices`, `faces_edges`, and `faces_faces` from this frame.
    pub fn extract_topology(self)
        -> Result<(Self, Option<TiVec<Vertex, Vec<Vertex>>>, Option<TiVec<Vertex, Vec<Edge>>>, Option<TiVec<Vertex, Vec<Option<Face>>>>,
                         Option<TiVec<Edge, [Vertex; 2]>>,                                     Option<TiVec<Edge, Vec<Option<Face>>>>,
                         Option<TiVec<Face, Vec<Vertex>>>,   Option<TiVec<Face, Vec<Edge>>>,   Option<TiVec<Face, Vec<Option<Face>>>>), String>
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
    use typed_index_collections::{ti_vec, TiVec};

    use crate::{filter::{edges_vertices_incident, other_vertex}, ser_de::{SerDeFrame, SerDeRat}, topology::TopologyError};
    use crate::fold::{Vertex as V, Edge as E, Face as F, HalfEdge as H};

    #[test]
    fn test_num_vertices_vertices_coords_f64() {
        let frame = SerDeFrame {
            vertices_coords_f64: Some(ti_vec![
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
            vertices_coords_exact: Some(ti_vec![
                vec![vec![SerDeRat(Rat::ONE)], vec![SerDeRat(Rat::TWO)]],
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((1, "vertices_exact:coords".to_owned()))));
    }

    #[test]
    fn test_num_vertices_vertices_vertices() {
        let frame = SerDeFrame {
            vertices_vertices: Some(ti_vec![
                vec![V(1), V(4), V(5)],
                vec![V(3), V(4), V(1), V(2)]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((2, "vertices_vertices".to_owned()))));
    }

    #[test]
    fn test_num_vertices_vertices_edges() {
        let frame = SerDeFrame {
            vertices_edges: Some(ti_vec![
                vec![E(1), E(4), E(5)],
                vec![E(3), E(4), E(1), E(2)]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((2, "vertices_edges".to_owned()))));
    }

    #[test]
    fn test_num_vertices_vertices_edges_mismatched() {
        let frame = SerDeFrame {
            vertices_vertices: Some(ti_vec![
                vec![V(1), V(2)],
            ]),
            vertices_edges: Some(ti_vec![
                vec![E(1), E(4), E(5)],
                vec![E(3), E(4), E(1), E(2)]
            ]),
            ..Default::default()
        };
        assert!(frame.num_vertices(&IndexMap::new()).is_err());
    }

    #[test]
    fn test_num_vertices_vertices_faces() {
        let frame = SerDeFrame {
            vertices_faces: Some(ti_vec![
                vec![Some(F(1)), None, Some(F(3))],
                vec![Some(F(0)), Some(F(3)), Some(F(4))],
                vec![None, Some(F(2)), Some(F(5))],
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
                ("test:field".to_owned(), ti_vec![
                    json!(3),
                    json!(5),
                    json!(5),
                ])
            ].into_iter().collect::<IndexMap<_, _>>()), Ok(Some((3, "vertices_test:field".to_owned()))));
    }

    #[test]
    fn test_num_vertices_edges_vertices() {
        let frame = SerDeFrame {
            edges_vertices: Some(ti_vec![
                [V(3), V(4)],
                [V(1), V(4)],
                [V(5), V(2)]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((6, "edges_vertices".to_owned()))));
    }

    #[test]
    fn test_num_vertices_faces_vertices() {
        let frame = SerDeFrame {
            faces_vertices: Some(ti_vec![
                vec![V(0), V(3), V(4)],
                vec![V(2), V(1), V(4)],
                vec![V(3), V(5), V(2)]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((6, "faces_vertices".to_owned()))));
    }

    #[test]
    fn test_num_vertices_edges_vertices_faces_vertices() {
        let mut frame = SerDeFrame {
            edges_vertices: Some(ti_vec![
                [V(1), V(4)],
                [V(2), V(3)],
            ]),
            faces_vertices: Some(ti_vec![
                vec![V(2), V(7), V(4)],
                vec![V(3), V(4), V(2)]
            ]),
            ..Default::default()
        };
        assert_eq!(frame.num_vertices(&IndexMap::new()), Ok(Some((8, "faces_vertices".to_owned()))));
        frame.edges_vertices = Some(ti_vec![
            [V(10), V(2)],
            [V(3), V(5)]
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
    fn test_try_ev_to_vh() {
        let vertices_edges = super::try_ev_to_vh(&ti_vec![], 0);
        assert_eq!(vertices_edges, Ok(ti_vec![]));

        // Star
        let mut vertices_edges = super::try_ev_to_vh(&ti_vec![
            [V(0), V(1)],
            [V(0), V(2)],
            [V(0), V(3)],
        ], 4);
        vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(vertices_edges, Ok(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)]
        ]));

        // Slash
        let mut vertices_edges = super::try_ev_to_vh(&ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(1), V(3)],
        ], 4);
        vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(vertices_edges, Ok(ti_vec![
            vec![H(0), H(7)],
            vec![H(1), H(2), H(8)],
            vec![H(3), H(4)],
            vec![H(5), H(6), H(9)],
        ]));

        // Doubled edge
        let mut vertices_edges = super::try_ev_to_vh(&ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(1), V(2)],
        ], 3);
        vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(vertices_edges, Ok(ti_vec![
            vec![H(0)],
            vec![H(1), H(2), H(4)],
            vec![H(3), H(5)],
        ]));

        // Isolated vertex
        let mut vertices_edges = super::try_ev_to_vh(&ti_vec![
            [V(0), V(1)],
            [V(0), V(2)],
        ], 4);
        vertices_edges.as_mut().unwrap().iter_mut().for_each(|v| v.sort());
        assert_eq!(vertices_edges, Ok(ti_vec![
            vec![H(0), H(2)],
            vec![H(1)],
            vec![H(3)],
            vec![],
        ]));

        // Self-loop
        let vertices_edges = super::try_ev_to_vh(&ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(2)],
        ], 3);
        assert!(vertices_edges.is_err());
    }

    fn canonicalize_vh_ev(mut vh_ev: (TiVec<V, Vec<H>>, TiVec<E, [V; 2]>)) -> (TiVec<V, Vec<H>>, TiVec<E, [V; 2]>) {
        for (e, vs) in vh_ev.1.iter_mut_enumerated() {
            if vs[0] > vs[1] {
                vs.reverse();
                vh_ev.0[vs[0]].iter_mut().find(|h| h.edge() == e).unwrap().flip();
                vh_ev.0[vs[1]].iter_mut().find(|h| h.edge() == e).unwrap().flip();
            }
        }
        vh_ev
    }

    #[test]
    fn test_try_ve_to_vh_ev() {
        use TopologyError::*;
        let edges_vertices = super::try_ve_to_vh_ev(&ti_vec![]);
        assert_eq!(edges_vertices, Ok((ti_vec![], ti_vec![])));

        // Star
        let edges_vertices = super::try_ve_to_vh_ev(&ti_vec![
            vec![E(0), E(1), E(2)],
            vec![E(0)],
            vec![E(1)],
            vec![E(2)]
        ]).map(canonicalize_vh_ev);
        assert_eq!(edges_vertices, Ok((ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)]
        ], ti_vec![
            [V(0), V(1)],
            [V(0), V(2)],
            [V(0), V(3)],
        ])).map(canonicalize_vh_ev));

        // Slash
        let edges_vertices = super::try_ve_to_vh_ev(&ti_vec![
            vec![E(0), E(3)],
            vec![E(0), E(1), E(4)],
            vec![E(1), E(2)],
            vec![E(2), E(3), E(4)],
        ]).map(canonicalize_vh_ev);
        assert_eq!(edges_vertices, Ok((ti_vec![
            vec![H(0), H(6)],
            vec![H(1), H(2), H(8)],
            vec![H(3), H(4)],
            vec![H(5), H(7), H(9)],
        ], ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(0), V(3)],
            [V(1), V(3)],
        ])).map(canonicalize_vh_ev));

        // Doubled edge
        let edges_vertices = super::try_ve_to_vh_ev(&ti_vec![
            vec![E(0)],
            vec![E(0), E(1), E(2)],
            vec![E(1), E(2)],
        ]).map(canonicalize_vh_ev);
        assert_eq!(edges_vertices, Ok((ti_vec![
            vec![H(0)],
            vec![H(1), H(2), H(4)],
            vec![H(3), H(5)],
        ], ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(1), V(2)],
        ])).map(canonicalize_vh_ev));

        // Isolated vertex
        let edges_vertices = super::try_ve_to_vh_ev(&ti_vec![
            vec![E(0), E(1)],
            vec![E(0)],
            vec![E(1)],
            vec![],
        ]).map(canonicalize_vh_ev);
        assert_eq!(edges_vertices, Ok((ti_vec![
            vec![H(0), H(2)],
            vec![H(1)],
            vec![H(3)],
            vec![],
        ], ti_vec![
            [V(0), V(1)],
            [V(0), V(2)],
        ])).map(canonicalize_vh_ev));

        // Self-loop
        let edges_vertices = super::try_ve_to_vh_ev(&ti_vec![
            vec![E(0)],
            vec![E(0), E(2)],
            vec![E(2), E(1), E(1)],
        ]);
        assert_eq!(edges_vertices, Err(vec![EdgeIsSelfLoop { edge: E(1), vertex: V(2) }]));

        // Edge index skipped
        let edges_vertices = super::try_ve_to_vh_ev(&ti_vec![
            vec![E(0), E(1)],
            vec![E(1), E(3)],
            vec![E(3), E(0)],
        ]);
        assert_eq!(edges_vertices, Err(vec![EdgeDoesNotHave2Vertices { edge: E(2), vertices: vec![] }]));

        // Edge missing vertices
        let edges_vertices = super::try_ve_to_vh_ev(&ti_vec![
            vec![E(0), E(1)],
            vec![E(2)],
            vec![E(2), E(0)],
        ]);
        assert_eq!(edges_vertices, Err(vec![EdgeDoesNotHave2Vertices { edge: E(1), vertices: vec![V(0)] }]));

        // Edge has too many vertices
        let mut edges_vertices = super::try_ve_to_vh_ev(&ti_vec![
            vec![E(0), E(1)],
            vec![E(1), E(2)],
            vec![E(2), E(0)],
            vec![E(1)],
        ]);
        let _ =
            edges_vertices.as_mut().map_err(|e| if let EdgeDoesNotHave2Vertices{ vertices, ..} = &mut e[0] { vertices.sort() } );
        assert_eq!(edges_vertices, Err(vec![EdgeDoesNotHave2Vertices { edge: E(1), vertices: vec![V(0), V(1), V(3)] }]));
    }

    #[test]
    fn test_try_fe_ev_to_fh() {
        use TopologyError::*;
        let faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![], &ti_vec![]);
        assert_eq!(faces_half_edges, Ok(ti_vec![]));

        // Single triangle
        let faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![
            vec![E(0), E(2), E(1)],
        ], &ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(0)],
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(1), H(5), H(3)]
        ]));

        // Doubly-covered square
        let faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![
            vec![E(0), E(1), E(2), E(3)],
            vec![E(3), E(2), E(1), E(0)]
        ], &ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)]
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ]));

        // Square cross
        let faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![
            vec![E(0), E(5), E(4)],
            vec![E(1), E(6), E(5)],
            vec![E(2), E(7), E(6)],
            vec![E(3), E(4), E(7)]
        ], &ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(4)],
            [V(1), V(4)],
            [V(2), V(4)],
            [V(3), V(4)]
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ]));

        // Slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![
            vec![E(0), E(1), E(2), E(3), E(8), E(7), E(6), E(5), E(4), E(8)],
            vec![E(4), E(5), E(6), E(7)],
        ], &ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(4), V(5)],
            [V(5), V(6)],
            [V(6), V(7)],
            [V(7), V(4)],
            [V(0), V(4)],
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
            vec![H(8), H(10), H(12), H(14)],
        ]));

        // Partially slitted face
        let faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![
            vec![E(4), E(4), E(0), E(1), E(2), E(3)],
        ], &ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(4)],
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(8), H(9), H(0), H(2), H(4), H(6)],
        ]));
        
        // Doubled-up edge
        let mut faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![
            vec![E(0), E(1), E(4)],
            vec![E(2), E(3), E(5)],
            vec![E(4), E(5)],
        ], &ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(2)],
            [V(0), V(2)],
        ]);
        let _ = faces_half_edges.as_mut().map(|fh| {
            if fh[F(2)][0].flip_bit() {
                fh[F(2)][0].flip();
                fh[F(2)][1].flip();
            }
        });
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(0), H(2), H(9)],
            vec![H(4), H(6), H(10)],
            vec![H(8), H(11)], // not unique!
            // vec![H(9), H(10)] This is also possible
        ]));

        // Face with 0 vertices. Not allowed.
        let faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![
            vec![E(0), E(1), E(2), E(3)],
            vec![],
        ], &ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
        ]);
        assert_eq!(faces_half_edges, Err(vec![
            FaceDoesNotHaveAtLeast2Edges { face: F(1), edges: vec![] }
        ]));

        // Open face. Not allowed.
        let faces_half_edges = super::try_fe_ev_to_fh(&ti_vec![
            vec![E(0), E(1), E(2)],
        ], &ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
        ]);
        assert_eq!(std::mem::discriminant(&faces_half_edges.unwrap_err()[0]),
            std::mem::discriminant(&ConsectiveEdgesOfFaceDoNotShareVertex { face: F(0), edges: [E(0); 2], vertices: [[V(0), V(0)]; 2] }),
        );
    }

    #[test]
    fn test_try_fv_vh_to_fh() {
        use TopologyError::*;
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![], &ti_vec![]);
        assert_eq!(faces_half_edges, Ok(ti_vec![]));

        // Single triangle
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![
            vec![V(1), V(0), V(2)],
        ], &ti_vec![
            vec![H(0), H(5)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(1), H(5), H(3)]
        ]));

        // Doubly-covered square
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![
            vec![V(0), V(1), V(2), V(3)],
            vec![V(0), V(3), V(2), V(1)],
        ], &ti_vec![
            vec![H(0), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ]));

        // Square cross
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![
            vec![V(0), V(1), V(4)],
            vec![V(1), V(2), V(4)],
            vec![V(2), V(3), V(4)],
            vec![V(3), V(0), V(4)],
        ], &ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(11), H(13), H(15)],
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ]));

        // Slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![
            vec![V(0), V(1), V(2), V(3), V(0), V(4), V(7), V(6), V(5), V(4)],
            vec![V(4), V(5), V(6), V(7)],
        ], &ti_vec![
            vec![H(0), H(7), H(16)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15), H(17)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
            vec![H(8), H(10), H(12), H(14)],
        ]));

        // Partially slitted face
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![
            vec![V(0), V(4), V(0), V(1), V(2), V(3)],
        ], &ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(9)],
        ]);
        assert_eq!(faces_half_edges, Ok(ti_vec![
            vec![H(8), H(9), H(0), H(2), H(4), H(6)],
        ]));
        
        // Doubled-up edge. Not allowed
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![
            vec![V(0), V(1), V(2)],
            vec![V(2), V(3), V(0)],
            vec![V(0), V(2)],
        ], &ti_vec![
            vec![H(0), H(7), H(8), H(10)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9), H(11)],
            vec![H(5), H(6)],
        ]);
        assert_eq!(std::mem::discriminant(&faces_half_edges.unwrap_err()[0]),
            std::mem::discriminant(&AmbiguousDoubleEdge { face: F(0), edges: vec![], vertices: [V(0), V(0)] }));

        // Face with 0 vertices. Not allowed.
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![
            vec![V(0), V(1), V(2), V(3)],
            vec![],
        ], &ti_vec![
            vec![H(0), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
        ]);
        assert_eq!(faces_half_edges, Err(vec![
            FaceDoesNotHaveAtLeast2Vertices { face: F(1), vertices: vec![] }
        ]));

        // Open face. Not allowed.
        let faces_half_edges = super::try_fv_vh_to_fh(&ti_vec![
            vec![V(0), V(1), V(2)],
        ], &ti_vec![
            vec![H(0), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
        ]);
        assert_eq!(std::mem::discriminant(&faces_half_edges.unwrap_err()[0]),
            std::mem::discriminant(&ConsectiveVerticesOfFaceDoNotShareEdge { face: F(0), vertices: [V(0), V(0)] }),
        );
    }

    fn canonicalize_ev_fh(mut ev_fh: (TiVec<E, [V; 2]>, TiVec<F, Vec<H>>)) -> (TiVec<E, [V; 2]>, TiVec<F, Vec<H>>) {
        for (e, vs) in ev_fh.0.iter_mut_enumerated() {
            if vs[0] > vs[1] {
                vs.reverse();
                ev_fh.1.iter_mut().flatten().filter(|h| h.edge() == e).for_each(|h| h.flip());
            }
        }
        ev_fh
    }

    #[test]
    fn test_try_fv_fe_to_ev_fh() {
        use TopologyError::*;
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![], &ti_vec![], 0);
        assert_eq!(ev_fh, Ok((ti_vec![], ti_vec![])));

        // Single triangle
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(1), V(0), V(2)],
        ], &ti_vec![
            vec![E(0), E(2), E(1)],
        ], 3).map(canonicalize_ev_fh);
        assert_eq!(ev_fh, Ok((ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(0)],
        ], ti_vec![
            vec![H(1), H(5), H(3)]
        ])).map(canonicalize_ev_fh));

        // Doubly-covered square
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(0), V(1), V(2), V(3)],
            vec![V(0), V(3), V(2), V(1)],
        ], &ti_vec![
            vec![E(0), E(1), E(2), E(3)],
            vec![E(3), E(2), E(1), E(0)]
        ], 4).map(canonicalize_ev_fh);
        assert_eq!(ev_fh, Ok((ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)]
        ], ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ])).map(canonicalize_ev_fh));

        // Square cross
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(0), V(1), V(4)],
            vec![V(1), V(2), V(4)],
            vec![V(2), V(3), V(4)],
            vec![V(3), V(0), V(4)],
        ], &ti_vec![
            vec![E(0), E(5), E(4)],
            vec![E(1), E(6), E(5)],
            vec![E(2), E(7), E(6)],
            vec![E(3), E(4), E(7)]
        ], 8).map(canonicalize_ev_fh);
        assert_eq!(ev_fh, Ok((ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(4)],
            [V(1), V(4)],
            [V(2), V(4)],
            [V(3), V(4)]
        ], ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ])).map(canonicalize_ev_fh));

        // Slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(0), V(1), V(2), V(3), V(0), V(4), V(7), V(6), V(5), V(4)],
            vec![V(4), V(5), V(6), V(7)],
        ], &ti_vec![
            vec![E(0), E(1), E(2), E(3), E(8), E(7), E(6), E(5), E(4), E(8)],
            vec![E(4), E(5), E(6), E(7)],
        ], 9).map(canonicalize_ev_fh);
        assert_eq!(ev_fh, Ok((ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(4), V(5)],
            [V(5), V(6)],
            [V(6), V(7)],
            [V(7), V(4)],
            [V(0), V(4)],
        ], ti_vec![
            vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
            vec![H(8), H(10), H(12), H(14)],
        ])).map(canonicalize_ev_fh));

        // Partially slitted face
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(0), V(4), V(0), V(1), V(2), V(3)],
        ], &ti_vec![
            vec![E(4), E(4), E(0), E(1), E(2), E(3)],
        ], 5).map(canonicalize_ev_fh);
        assert_eq!(ev_fh, Ok((ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(4)],
        ], ti_vec![
            vec![H(8), H(9), H(0), H(2), H(4), H(6)],
        ])).map(canonicalize_ev_fh));
        
        // Doubled-up edge
        let mut ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(0), V(1), V(2)],
            vec![V(2), V(3), V(0)],
            vec![V(0), V(2)],
        ], &ti_vec![
            vec![E(0), E(1), E(4)],
            vec![E(2), E(3), E(5)],
            vec![E(4), E(5)],
        ], 6).map(canonicalize_ev_fh);
        assert_eq!(ev_fh, Ok((ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(2)],
            [V(0), V(2)],
        ], ti_vec![
            vec![H(0), H(2), H(9)],
            vec![H(4), H(6), H(10)],
            vec![H(8), H(11)],
        ])).map(canonicalize_ev_fh));

        // Face with 0 vertices. Not allowed.
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(0), V(1), V(2), V(3)],
            vec![],
        ], &ti_vec![
            vec![E(0), E(1), E(2), E(3)],
            vec![],
        ], 4);
        assert_eq!(ev_fh, Err(vec![
            FaceDoesNotHaveAtLeast2Vertices { face: F(1), vertices: vec![] }
        ]));

        // Face with 1 vertex. Not allowed.
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(0), V(1), V(2), V(3)],
            vec![V(5)],
        ], &ti_vec![
            vec![E(0), E(1), E(2), E(3)],
            vec![E(4)],
        ], 5);
        assert_eq!(ev_fh, Err(vec![
            FaceDoesNotHaveAtLeast2Vertices { face: F(1), vertices: vec![V(5)] }
        ]));

        // Faces and edges don't have the same counts.
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(1), V(0), V(2)],
        ], &ti_vec![
            vec![E(0), E(2)],
        ], 3);
        assert_eq!(ev_fh.unwrap_err()[0], 
            FaceVertexEdgeCountMismatch { face: F(0), vertices: vec![V(1), V(0), V(2)], edges: vec![E(0), E(2)] }
        );

        // Edge definition conflict
        let ev_fh = super::try_fv_fe_to_ev_fh(&ti_vec![
            vec![V(1), V(0), V(2)],
            vec![V(3), V(4)],
        ], &ti_vec![
            vec![E(0), E(2), E(1)],
            vec![E(0), E(5)],
        ], 6);
        assert_eq!(std::mem::discriminant(&ev_fh.unwrap_err()[0]),
            std::mem::discriminant(&EdgeVerticesMismatch { face: F(0), edge: E(0), vertices: [[V(0); 2]; 2] }),
        );
    }
}