use std::fmt::Display;

use indexmap::{IndexSet};
use typed_index_collections::{ti_vec, TiVec};

use crate::fold::{AtFaceCorner, EdgesFaceCornersEx, Edge, Face, FaceCorner, Frame, FrameAttribute, HalfEdge, Vertex};

/// An error that happens when checking if a frame is a manifold
#[derive(Clone, Debug)]
#[cfg_attr(test, derive(PartialEq, Eq, PartialOrd, Ord))] // Really no point in this derivation outside of tests
pub enum ManifoldError {
    /// An edge is adjacent to at least 3 faces, which is too many.
    EdgeAdjacentToTooManyFaceCorners { edge: Edge, face_corners: [Vec<FaceCorner>; 2] },
    /// A vertex's face corners cannot be sequenced in a loop or some fans
    VertexFacesDoNotFormFansOrLoop { vertex: Vertex, face_corners: Vec<FaceCorner> }
}

impl Display for ManifoldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "manifold error: ")?;
        match self {
            Self::EdgeAdjacentToTooManyFaceCorners { edge, face_corners } => {
                write!(f, "edge {edge} is adjacent to 3 or more face corners ({face_corners:?})")
            }
            Self::VertexFacesDoNotFormFansOrLoop { vertex, face_corners } => {
                write!(f, "vertex {vertex}'s face corners ({face_corners:?}) do not form a fan or loop")
            }
        }
    }
}

/// An error that happens when checking if a frame is an orientable manifold
#[derive(Clone, Debug)]
#[cfg_attr(test, derive(PartialEq, Eq, PartialOrd, Ord))] // Really no point in this derivation outside of tests
pub enum OrientableError {
    NotAManifold(ManifoldError),
    /// The following faces cannot be oriented consistently.
    CannotOrientFacesConsistently { faces: Vec<Face> },
    /// A vertex's face corners cannot be sequenced in a loop or some fans
    VertexFacesDoNotFormFansOrLoop { vertex: Vertex, face_corners: Vec<FaceCorner> }
}

impl Display for OrientableError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "manifold error: ")?;
        match self {
            Self::NotAManifold(e) => {
                write!(f, "{e}")
            }
            Self::CannotOrientFacesConsistently { faces } => {
                write!(f, "the faces {faces:?} cannot be oriented consistently")
            }
            Self::VertexFacesDoNotFormFansOrLoop { vertex, face_corners } => {
                write!(f, "vertex {vertex}'s face corners ({face_corners:?}) do not form a fan or loop")
            }
        }
    }
}

impl Frame {
    /// Gets a partition of fans of the vertex.
    /// 
    /// A *fan* is a sequence of half-edges where the corresponding edges of consecutive half-edges share a face.
    /// If the edges of the first and last half-edges also share a face, then the fan is a *loop*.
    /// 
    /// The returned list is a vector of tuples of (fan, is_loop).
    /// Returns `None` if any half-edge's edge is adjacent to more than 2 faces,
    /// because then the concept of a fan partition breaks down.
    pub fn vertex_fans(&self, vertex: Vertex) -> Option<Vec<(Vec<HalfEdge>, bool)>> {
        // Wow this is complicated.
        let vh = if let Some(vh) = self.vertices_half_edges.as_ref() { vh } else {
            return Some(vec![]);
        };
        let ec = if let Some(ec) = self.edges_face_corners.as_ref() { ec } else {
            return Some(vh[vertex].iter().map(|&h| (vec![h], false)).collect());
        };
        let fh = self.faces_half_edges.as_ref().unwrap(); // exists if `edges_face_corners` exists

        if vh[vertex].iter().any(|&h| !self.is_edge_manifold(h.edge())) {
            return None;
        }

        let mut unvisited_h = vh[vertex].iter().copied().collect::<IndexSet<_>>();
        let mut fans = vec![];

        while !unvisited_h.is_empty() {
            // Find a start: an edge adjacent to exactly 1 face.
            // Otherwise, just take an edge adjacent to 2 faces; otherwise skip
            let start = unvisited_h.iter().find(|&&h|
                ec.at(h).len() + ec.at(h.flipped()).len() == 1).copied();
            let start = start.or_else(|| unvisited_h.iter().find(|&&h|
                ec.at(h).len() + ec.at(h.flipped()).len() == 2).copied());
            let start_half_edge = if let Some(s) = start { s } else {
                fans.extend(unvisited_h.into_iter().map(|h| (vec![h], false)));
                return Some(fans);
            };

            // Build the fan
            let mut curr_half_edge = start_half_edge;
            let mut curr_face_corner = None;
            let mut fan = vec![];
            let mut is_loop = false;
            loop {
                fan.push(curr_half_edge);
                unvisited_h.swap_remove(&curr_half_edge);

                // Next half-edge
                let corner = ec.pair_at(curr_half_edge).iter().enumerate()
                    .flat_map(|(flip, cs)|
                        cs.iter().map(move |&c| (flip != 0, if flip != 0 { c.next(fh) } else { c })))
                    .find(|&(_, c)| curr_face_corner.map(|c2| c != c2).unwrap_or(true));
                let (flip, corner) = if let Some(c) = corner { c } else { break };
                curr_half_edge = if flip {
                    fh.at(corner)
                } else {
                    fh.at(corner.prev(fh)).flipped()
                };

                curr_face_corner = Some(corner);

                if curr_half_edge == start_half_edge {
                    is_loop = true;
                    break;
                }
            }

            fans.push((fan, is_loop));
        }
        Some(fans)
    }
    
    /// Gets a partition of fans of the vertex if the vertex is manifold.
    /// See `Frame::vertex_fans`
    pub fn vertex_fans_if_manifold(&self, vertex: Vertex) -> Option<Vec<(Vec<HalfEdge>, bool)>> {
        let fans = if let Some(f) = self.vertex_fans(vertex) { f } else { return None };
        if fans.len() <= 1 || fans.iter().all(|(_, is_loop)| !is_loop) { Some(fans) } else { None }
    }

    /// Checks if a vertex is manifold,
    /// using the multiple open fans or one loop definition.
    pub fn is_vertex_manifold(&self, vertex: Vertex) -> bool {
        self.vertex_fans_if_manifold(vertex).is_some()
    }

    /// Checks if an edge is manifold.
    /// An edge is manifold if it is used by at most 2 face corners.
    /// (not faces, *face corners*, to account for the fact that sometimes
    /// an edge is used multiple times by the same face.)
    pub fn is_edge_manifold(&self, edge: Edge) -> bool {
        // *much* simpler than that vertex definition
        let ec = if let Some(ec) = self.edges_face_corners.as_ref() { ec } else { return true };
        ec[edge][0].len() + ec[edge][1].len() <= 2
    }

    /// Checks if an edge is oriented.
    /// An edge is oriented if it is used by at most 2 face corners,
    /// and no two face corners wind it in the same direction.
    pub fn is_edge_oriented(&self, edge: Edge) -> bool {
        let ec = if let Some(ec) = self.edges_face_corners.as_ref() { ec } else { return true };
        ec[edge][0].len() <= 1 && ec[edge][1].len() <= 1
    }

    /// Tries to convert this frame into a manifold (possibly with boundary) representation.
    /// On success, every vertex will have its half-edges listed cyclically according to the fans/loop of its faces,
    /// and `FrameAttribute::Manifold` is added to the frame attributes. Nothing else changes.
    /// 
    /// The result includes for each vertex, the delimiter indices of its fans. If a vertex has a loop,
    /// that list is empty since there are no delimiters. Otherwise, the delimiter list includes 0 and the number of half-edges
    /// in the vertex.
    /// 
    /// In the case that a vertex is adjacent to multiple fans,
    /// the half-edges of each fan are consecutive, and the order of the fans (and orientation of each fan) is arbitrary
    /// 
    /// Requirements:
    /// * Each edge must be adjacent to at most 2 faces
    /// * TODO: Define the vertex fan condition, the below is wrong
    /// * ~~In each vertex, the faces must be sequenceable in such a way that adjacent faces share an edge.
    ///     (in other words, the faces form a fan or loop. Double-cones are not allowed.)~~
    /// 
    /// Note that a frame with no faces is always a manifold according to this definition.
    pub fn try_into_manifold(mut self) -> Result<(Self, TiVec<Vertex, Vec<usize>>), Vec<ManifoldError>> {
        self.add_attribute_unchecked(FrameAttribute::Manifold);
        self.remove_attribute(FrameAttribute::NonManifold);

        // Welp; now to check it.
        let (vh, ef) = if let (Some(vh), Some(ef)) =
            (self.vertices_half_edges.as_mut(), self.edges_face_corners.as_ref())
        {
            (vh, ef)
        } else {
            let fan_delimiters = if let Some(vh) = self.vertices_half_edges.as_ref() {
                vh.iter().map(|hs| (0..=hs.len()).collect::<Vec<_>>()).collect::<TiVec<_, _>>()
            } else {
                ti_vec![]
            };
            return Ok((self, fan_delimiters))
        };

        let mut errors = vec![];
        for (e, faces) in ef.iter_enumerated() {
            if faces[0].len() + faces[1].len() > 2 {
                errors.push(ManifoldError::EdgeAdjacentToTooManyFaceCorners { edge: e, face_corners: faces.clone() })
            }
        }

        let mut vertices_fan_boundaries = ti_vec![vec![]; vh.len()];

        for v in (0..vh.len()).map(Vertex) {
            let fans = if let Some(f) = self.vertex_fans_if_manifold(v) { f } else {
                let vh = self.vertices_half_edges.as_ref().unwrap();
                errors.push(ManifoldError::VertexFacesDoNotFormFansOrLoop {
                    vertex: v,
                    face_corners: vh[v].iter().map(|h| *h)
                        .flat_map(|h| ef.at(h).iter()).copied().collect::<Vec<_>>()
                });
                continue;
            };
            let vh = self.vertices_half_edges.as_mut().unwrap();

            vh[v].clear();
            let is_loop = fans.len() == 1 && fans[0].1;
            for (fan, _) in fans {
                vertices_fan_boundaries[v].push(vh[v].len());
                vh[v].extend(fan);
            }
            vertices_fan_boundaries[v].push(vh[v].len());
            if is_loop {
                vertices_fan_boundaries[v].clear();
            }
        }

        if errors.is_empty() { Ok((self, vertices_fan_boundaries)) } else { Err(errors) }
    }

    /// Tries to convert this frame into an orientable manifold (possibly with boundary) representation.
    /// On success:
    /// * every vertex will have its half-edges listed counterclockwise according to the fan/loop of its faces
    ///     (however, the manifold could be inside-out, turning counterclockwise into clockwise)
    /// * some faces may get their orientation flipped for consistency, affecting `faces_half_edges` and `edges_face_corners`.
    ///     (specifically, the list of half-edges in a face could get reversed. It will not be permuted in a different way,
    ///     not even a cyclic permutation.)
    /// * each connected-component-by-face will contain a face that has *not* been flipped. In particular, if a connected-component-by-face
    ///     was already oriented, no faces in it will be flipped.
    /// * `FrameAttribute::Manifold` and `FrameAttribute::Orientable` get added to the attribute list
    /// 
    /// The result includes for each vertex, the delimiter indices of its fans. If a vertex has a loop,
    /// that list is empty since there are no delimiters. Otherwise, the delimiter list includes 0 and the number of half-edges
    /// in the vertex.
    /// 
    /// In the case that a vertex is adjacent to some edges with 0 faces attached to them,
    /// the corresponding half-edges go to the end of the list in their previous order.
    /// 
    /// Requirements:
    /// * The faces must be reorientable in such a way that each *half-edge* is adjacent to at most 1 face
    ///     (in other words, such that each *edge* is adjacent to at most 1 face winding it forward and at most 1 face winding it backward)
    /// * In each vertex, the faces must be sequenceable in such a way that adjacent faces share an edge.
    ///     (in other words, the faces form a fan or loop. Double-cones are not allowed.)
    /// 
    /// Note that a frame with no faces is always a manifold according to this definition.
    pub fn try_into_orientable(mut self) -> Result<(Self, TiVec<Vertex, Vec<usize>>), Vec<OrientableError>> {
        //TODO: Test. This is extremely annoying to test, unfortunately.
        let (new_self, fan_delimiters) = self.try_into_manifold()
            .map_err(|e| e.into_iter().map(OrientableError::NotAManifold).collect::<Vec<_>>())?;
        self = new_self;

        self.add_attribute_unchecked(FrameAttribute::Orientable);
        self.remove_attribute(FrameAttribute::NonOrientable);

        // Welp; now to check it.
        let (vh, _ev, ef, fh) = if let (Some(vh), Some(ev), Some(ef), Some(fh)) =
            (self.vertices_half_edges.as_mut(), self.edges_vertices.as_ref(), self.edges_face_corners.as_mut(), self.faces_half_edges.as_mut())
        {
            (vh, ev, ef, fh)
        } else { return Ok((self, fan_delimiters)) };

        // Attempt to orient all faces
        let mut unvisited_faces = fh.iter_enumerated().map(|(f, _)| f).collect::<IndexSet<_>>();

        // Beware of multiple connected components
        while let Some(&unvisited) = unvisited_faces.first() {
            let mut faces_to_visit = std::iter::once(unvisited).collect::<IndexSet<_>>();
            let mut this_component = vec![];

            while let Some(f) = faces_to_visit.pop() {
                if !unvisited_faces.swap_remove(&f) { continue; }
                this_component.push(f);
                
                if fh[f].iter().enumerate().any(|(i, &h)|
                    ef.at(h).iter().any(|&c| c != FaceCorner(f, i) && !unvisited_faces.contains(&c.face())))
                {
                    // Bad edge, flip the whole face!
                    let face_len = fh[f].len();
                    fh[f].iter_mut().enumerate().for_each(|(i, h)| {
                        let pos = ef.at(*h).iter().position(|&c| c == FaceCorner(f, i)).unwrap();
                        ef.at_mut(*h).swap_remove(pos);
                        ef.at_mut(h.flipped()).push(FaceCorner(f, face_len - i - 1));
                        h.flip();
                    });
                    fh[f].reverse();

                    if fh[f].iter().enumerate().any(|(i, &h)|
                        ef.at(h).iter().any(|&c| c != FaceCorner(f, i) && !unvisited_faces.contains(&c.face())))
                    {
                        return Err(vec![OrientableError::CannotOrientFacesConsistently { faces: this_component }]);
                    }
                }

                faces_to_visit.extend(fh[f].iter().flat_map(|&h| ef[h.edge()].iter().flatten()).map(|c| c.face()));
            }
        }

        // Order the fans counterclockwise. Note that we already know it's a manifold.
        for (half_edges, fans) in vh.iter_mut().zip(fan_delimiters.iter()) {
            if fans.is_empty() {
                // It's a loop. Make sure it's oriented counterclockwise
                if half_edges.len() >= 2 && fh.at(ef.at(half_edges[0])[0].prev(&fh)).flipped() != half_edges[1] {
                    half_edges.reverse();
                }
            } else {
                for (&fan_start, &fan_end) in fans.iter().zip(fans.iter().skip(1)) {
                    if !ef.at(half_edges[fan_start].flipped()).is_empty() {
                        half_edges[fan_start..fan_end].reverse();
                    }
                }
            }
        }

        Ok((self, fan_delimiters))
    }
}

#[cfg(test)]
mod test {
    use indexmap::indexset;
    use petgraph::{Graph};
    use typed_index_collections::{ti_vec, TiSlice, TiVec};

    use crate::fold::{Frame, FrameAttribute};
    use crate::fold::{Vertex as V, Edge as E, Face as F, HalfEdge as H};
    use crate::manifold::{ManifoldError, OrientableError};
    use crate::test_utils::assert_ec_fh_consistent;

    fn manifold_vh_graph(vh: &TiSlice<V, Vec<H>>, fan_delimiters: &TiSlice<V, Vec<usize>>) -> Graph<(V, H), ()> {
        let mut graph = Graph::new();
        let vertex_nodes = vh.iter_enumerated()
            .map(|(v, hs)| hs.iter()
                .map(|&h| graph.add_node((v, h)))
                .collect::<Vec<_>>())
            .collect::<TiVec<_, _>>();

        for (v, half_edges) in vh.iter_enumerated() {
            for (i, _) in half_edges.iter().enumerate() {
                if !fan_delimiters[v].contains(&(i + 1)) {
                    graph.add_edge(vertex_nodes[v][i], vertex_nodes[v][(i + 1) % vertex_nodes[v].len()], ());
                    graph.add_edge(vertex_nodes[v][(i + 1) % vertex_nodes[v].len()], vertex_nodes[v][i], ());
                }
            }
        }

        graph
    }

    #[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
    enum OrientableVhFhNodeLabel {
        Vertex(V),
        HalfEdge(H),
        Face(F),
    }

    fn orientable_vh_fh_graph(vh: &TiSlice<V, Vec<H>>, fh: &TiSlice<F, Vec<H>>, fan_delimiters: &TiSlice<V, Vec<usize>>)
        -> Graph<OrientableVhFhNodeLabel, ()>
    {
        use OrientableVhFhNodeLabel::*;

        let mut graph = Graph::new();

        // Vertex nodes and face nodes need duplication to handle the fact that you can turn each component inside-out
        let mut vertex_nodes = || vh.iter_enumerated()
            .map(|(v, hs)| hs.iter()
                .map(|_| graph.add_node(Vertex(v)))
                .collect::<Vec<_>>())
            .collect::<TiVec<_, _>>();
        let vertex_nodes_0 = vertex_nodes();
        let vertex_nodes_1 = vertex_nodes();

        let num_half_edges = vh.iter().flatten().copied().max().map(|h| h.0 + 1).unwrap_or(0);
        let mut half_edge_nodes = || (0..num_half_edges).map(|h| graph.add_node(HalfEdge(H(h)))).collect::<TiVec<_, _>>();
        let half_edge_nodes_0 = half_edge_nodes();
        let half_edge_nodes_1 = half_edge_nodes();

        let mut face_nodes = || fh.iter_enumerated()
            .map(|(f, hs)| hs.iter()
                .map(|_| graph.add_node(Face(f)))
                .collect::<Vec<_>>())
            .collect::<TiVec<_, _>>();
        let face_nodes_0 = face_nodes();
        let face_nodes_1 = face_nodes();

        // Connect vertices
        for (v, half_edges) in vh.iter_enumerated() {
            for (i, &h) in half_edges.iter().enumerate() {
                // to half edge
                graph.add_edge(vertex_nodes_0[v][i], half_edge_nodes_0[h], ());
                graph.add_edge(vertex_nodes_1[v][i], half_edge_nodes_1[h], ());
                if !fan_delimiters[v].contains(&(i + 1)) {
                    // in a chain to make a fan
                    graph.add_edge(vertex_nodes_0[v][i], vertex_nodes_0[v][(i + 1) % vertex_nodes_0[v].len()], ());
                    graph.add_edge(vertex_nodes_1[v][(i + 1) % vertex_nodes_1[v].len()], vertex_nodes_1[v][i], ()); // note the difference in copies
                }
            }
        }

        // Connect faces
        for (f, half_edges) in fh.iter_enumerated() {
            for (i, &h) in half_edges.iter().enumerate() {
                // to half edge
                graph.add_edge(face_nodes_0[f][i], half_edge_nodes_0[h], ());
                graph.add_edge(face_nodes_1[f][i], half_edge_nodes_1[h.flipped()], ());
                if i + 1 < half_edges.len() {
                    // in a chain to make a fan
                    graph.add_edge(face_nodes_0[f][i + 0], face_nodes_0[f][i + 1], ());
                    graph.add_edge(face_nodes_1[f][i + 1], face_nodes_1[f][i + 0], ()); // note the difference in copies
                }
            }
        }

        graph
    }

    fn assert_manifold_vh_isomorphic(a: &TiSlice<V, Vec<H>>, b: &TiSlice<V, Vec<H>>, fan_delimiters: &TiSlice<V, Vec<usize>>) {
        let iso = vf2::isomorphisms(
            &manifold_vh_graph(&a, &fan_delimiters), &manifold_vh_graph(&b, &fan_delimiters))
            .default_eq()
            .first();    

        assert!(iso.is_some(), "{a:?} is not isomorphic to {b:?}");
    }

    fn assert_manifold_vh_not_isomorphic(a: &TiSlice<V, Vec<H>>, b: &TiSlice<V, Vec<H>>, fan_delimiters: &TiSlice<V, Vec<usize>>) {
        let iso = vf2::isomorphisms(
            &manifold_vh_graph(&a, &fan_delimiters), &manifold_vh_graph(&b, &fan_delimiters))
            .default_eq()
            .first();    

        assert!(iso.is_none(), "{a:?} is isomorphic to {b:?}");
    }

    fn assert_orientable_vh_isomorphic(
        frame: &Frame,
        expected_vh: &TiSlice<V, Vec<H>>,
        expected_fh: &TiSlice<F, Vec<H>>,
        fan_delimiters: &TiSlice<V, Vec<usize>>
    ) {
        let default_vh = ti_vec![];
        let default_fh = ti_vec![];
        let default_ec = ti_vec![];
        let frame_vh = frame.vertices_half_edges.as_ref().unwrap_or(&default_vh);
        let frame_fh = frame.faces_half_edges.as_ref().unwrap_or(&default_fh);
        let frame_ec = frame.edges_face_corners.as_ref().unwrap_or(&default_ec);
        assert_ec_fh_consistent(frame_ec, frame_fh);

        let iso = vf2::isomorphisms(
            &orientable_vh_fh_graph(frame_vh, frame_fh, &fan_delimiters),
            &orientable_vh_fh_graph(&expected_vh, &expected_fh, &fan_delimiters)
        )
            .default_eq()
            .first();    

        assert!(iso.is_some(), "{frame_vh:?}\n{frame_fh:?}\nis not isomorphic to\n{expected_vh:?}\n{expected_fh:?}");
    }

    #[test]
    fn test_try_into_manifold() {
        use ManifoldError::*;

        // Empty frame
        let manifold = Frame {
            ..Default::default()
        }.try_into_manifold();
        let (manifold, _) = manifold.unwrap();
        assert_eq!(manifold.frame_attributes, indexset![FrameAttribute::Manifold]);
        assert_eq!(manifold.vertices_half_edges, None);

        // No faces. Just a 3-star.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ]), None).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ], &fans);

        // A simple triangle.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ], &fans);

        // A simple triangle, with some edges sticking out.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(6), H(10), H(8)],
            vec![H(7)],
            vec![H(9)],
            vec![H(11)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(6), H(10), H(8)],
            vec![H(7)],
            vec![H(9)],
            vec![H(11)],
        ], &fans);

        // A doubled square
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)]
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
        ], &fans);

        // A square cross. Note that the order of vertices in the apex matters a bit.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(13), H(11), H(15)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(7), H(8) , H(0)],
            vec![H(1), H(10), H(2)],
            vec![H(3), H(12), H(4)],
            vec![H(5), H(14), H(6)],
            vec![H(9), H(11), H(13), H(15)],
        ], &fans);
        assert_manifold_vh_not_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(7), H(8) , H(0)],
            vec![H(1), H(10), H(2)],
            vec![H(3), H(12), H(4)],
            vec![H(5), H(14), H(6)],
            vec![H(9), H(13), H(11), H(15)], // not allowed
        ], &fans);

        // A slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let manifold = Frame {
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
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(0), H(16), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15), H(17)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ], &fans);

        // A möbius face. Noticable, a face duplicates a half-edge.
        // Edges go 0->1->2->3->0 and 0->2.
        // Believe it or not, it's a manifold.
        // »---0       3--3»
        //     |"0   ."     
        //     4  "."       
        //     | 2" ".      
        // >---2"     "1--1>  
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9), H(7), H(5), H(9)],
            vec![H(0), H(2), H(4), H(6)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6)],
        ], &fans);

        // A möbius strip. This time, no face by itself violates orientability, but the whole does. Also, there's a boundary.
        // Edges go 0->1->2->3->0 and 0->2 and 3->1.
        // Believe it or not, it's a manifold.
        // »---0       3--3»
        //     |"0   ."|    
        //     4  "."  5    
        //     | 2" ". |    
        // >---2"     "1--1>  
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2), H(11)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6), H(10)],
        ]), Some(ti_vec![
            vec![H(0), H(11), H(5), H(9)],
            vec![H(6), H(8), H(3), H(11)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(7), H(8) , H(0)],
            vec![H(1), H(11), H(2)],
            vec![H(3), H(9) , H(4)],
            vec![H(5), H(10), H(6)],
        ], &fans);

        // A doubled-up edge
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8), H(10)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9), H(11)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9)],
            vec![H(4), H(6), H(10)],
            vec![H(8), H(11)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(0), H(8), H(10), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(9), H(11), H(4)],
            vec![H(5), H(6)],
        ], &fans);

        // A square cross, but two opposite faces are missing and there's subdivision happening
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(15)],
            vec![H(1), H(2), H(16)],
            vec![H(4), H(3)],
            vec![H(7), H(8)],
            vec![H(9), H(10), H(18)],
            vec![H(11), H(12)],
            vec![H(5), H(6), H(13), H(14), H(17), H(19)],
        ]), Some(ti_vec![
            vec![H(0), H(16), H(14)],
            vec![H(2), H(4), H(17)],
            vec![H(8), H(18), H(6)],
            vec![H(10), H(12), H(19)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(0), H(15)],
            vec![H(1), H(16), H(2)],
            vec![H(4), H(3)],
            vec![H(7), H(8)],
            vec![H(9), H(18), H(10)],
            vec![H(11), H(12)],
            vec![H(5), H(17), H(14), H(6), H(19), H(13)],
        ], &fans);

        // You know what's cooler than a square cross? *2* square crosses.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(13), H(11), H(15)],
            vec![H(16), H(23), H(24)],
            vec![H(17), H(18), H(26)],
            vec![H(19), H(20), H(28)],
            vec![H(21), H(22), H(30)],
            vec![H(25), H(29), H(27), H(31)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
            vec![H(16), H(26), H(25)],
            vec![H(18), H(28), H(27)],
            vec![H(20), H(30), H(29)],
            vec![H(22), H(24), H(31)],
        ])).try_into_manifold();
        let (manifold, fans) = manifold.unwrap();
        assert_manifold_vh_isomorphic(manifold.vertices_half_edges.as_ref().unwrap(), &ti_vec![
            vec![H(7), H(8) , H(0)],
            vec![H(1), H(10), H(2)],
            vec![H(3), H(12), H(4)],
            vec![H(5), H(14), H(6)],
            vec![H(9), H(11), H(13), H(15)],
            vec![H(16), H(24), H(23)],
            vec![H(17), H(26), H(18)],
            vec![H(19), H(28), H(20)],
            vec![H(21), H(30), H(22)],
            vec![H(25), H(27), H(29), H(31)],
        ], &fans);

        // Three triangles from a single edge (not a manifold)
        let not_a_manifold = Frame {
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
        ])).try_into_manifold();
        assert!(matches!(not_a_manifold.unwrap_err()[0], EdgeAdjacentToTooManyFaceCorners { edge: E(0), face_corners: _ }));

        // A square cross with an edge sticking out from the middle. Not a manifold
        let not_a_manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(13), H(11), H(15), H(16)],
            vec![H(17)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ])).try_into_manifold();
        assert!(matches!(not_a_manifold.unwrap_err()[0], VertexFacesDoNotFormFansOrLoop { vertex: V(4), face_corners: _ }));
    }

    #[test]
    fn test_try_into_orientable() {
        use OrientableError::*;

        // This is annoying to test; especially because you can turn connected components inside-out,
        // so just print results.

        // Empty frame
        let orientable = Frame {
            ..Default::default()
        }.try_into_orientable();
        let (orientable, _) = orientable.unwrap();
        assert_eq!(orientable.frame_attributes, indexset![FrameAttribute::Manifold, FrameAttribute::Orientable]);
        assert_eq!(orientable.vertices_half_edges, None);
 
        // No faces. Just a 3-star.
        let orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ]), None).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ], &ti_vec![], &fans);
        assert_eq!(orientable.faces_half_edges, None);
        assert_eq!(orientable.edges_face_corners, None);

        // A simple triangle.
        let orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ])).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ], &ti_vec![
            vec![H(5), H(3), H(1)]
        ], &fans);

        // A simple triangle, with some edges sticking out.
        let orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(6), H(10), H(8)],
            vec![H(7)],
            vec![H(9)],
            vec![H(11)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ])).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(6), H(10), H(8)],
            vec![H(7)],
            vec![H(9)],
            vec![H(11)],
        ], &ti_vec![
            vec![H(5), H(3), H(1)],
        ], &fans);

        // A doubled square
        let orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)]
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ])).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
        ], &ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ], &fans);

        // A square cross. Note that the order of vertices in the apex matters a bit.
        let orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(13), H(11), H(15)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(10), H(13), H(3)], // flip some faces to
            vec![H(12), H(15), H(5)], // make the algorithm do some work
            vec![H(6), H(8), H(15)],
        ])).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(0), H(8) , H(7)],
            vec![H(2), H(10), H(1)],
            vec![H(4), H(12), H(3)],
            vec![H(6), H(14), H(5)],
            vec![H(9), H(11), H(13), H(15)],
        ], &ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ], &fans);


        // A slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let orientable = Frame {
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
        ])).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(0), H(16), H(7)],
            vec![H(2), H(1)],
            vec![H(4), H(3)],
            vec![H(6), H(5)],
            vec![H(8), H(15), H(17)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ], &ti_vec![
            vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
            vec![H(8), H(10), H(12), H(14)],
        ], &fans);


        // A möbius face. Noticable, a face duplicates a half-edge.
        // Edges go 0->1->2->3->0 and 0->2.
        // Believe it or not, it's a orientable.
        // »---0       3--3»
        //     |"0   ."     
        //     4  "."       
        //     | 2" ".      
        // >---2"     "1--1>  
        let not_a_orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9), H(7), H(5), H(9)],
            vec![H(0), H(2), H(4), H(6)],
        ])).try_into_orientable();
        assert!(matches!(not_a_orientable.unwrap_err()[0], CannotOrientFacesConsistently { faces: _ }));

        // A möbius strip. This time, no face by itself violates orientability, but the whole does. Also, there's a boundary.
        // Edges go 0->1->2->3->0 and 0->2 and 3->1.
        // Believe it or not, it's a orientable.
        // »---0       3--3»
        //     |"0   ."|    
        //     4  "."  5    
        //     | 2" ". |    
        // >---2"     "1--1>  
        let not_a_orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2), H(11)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6), H(10)],
        ]), Some(ti_vec![
            vec![H(0), H(11), H(5), H(9)],
            vec![H(6), H(8), H(3), H(11)],
        ])).try_into_orientable();
        assert!(matches!(not_a_orientable.unwrap_err()[0], CannotOrientFacesConsistently { faces: _ }));

        // A doubled-up edge
        let orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8), H(10)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9), H(11)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9)],
            vec![H(4), H(6), H(10)],
            vec![H(8), H(11)],
        ])).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(7), H(10), H(8), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(9), H(11), H(4)],
            vec![H(5), H(6)],
        ], &ti_vec![
            vec![H(8), H(3), H(1)],
            vec![H(11), H(7), H(5)],
            vec![H(10), H(9)],
        ], &fans);

        // A square cross, but two opposite faces are missing and there's subdivision happening
        let orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(15)],
            vec![H(1), H(2), H(16)],
            vec![H(4), H(3)],
            vec![H(7), H(8)],
            vec![H(9), H(10), H(18)],
            vec![H(11), H(12)],
            vec![H(5), H(6), H(13), H(14), H(17), H(19)],
        ]), Some(ti_vec![
            vec![H(0), H(16), H(14)],
            vec![H(2), H(4), H(17)],
            vec![H(8), H(18), H(6)],
            vec![H(10), H(12), H(19)],
        ])).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(0), H(15)],
            vec![H(2), H(16), H(1)],
            vec![H(4), H(3)],
            vec![H(8), H(7)],
            vec![H(10), H(18), H(9)],
            vec![H(12), H(11)],
            vec![H(14), H(17), H(5), H(6), H(19), H(13)],
        ], &ti_vec![
            vec![H(0), H(16), H(14)],
            vec![H(2), H(4), H(17)],
            vec![H(8), H(18), H(6)],
            vec![H(10), H(12), H(19)],
        ], &fans);

        // You know what's cooler than a square cross? *2* square crosses.
        let orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(13), H(11), H(15)],
            vec![H(16), H(23), H(24)],
            vec![H(17), H(18), H(26)],
            vec![H(19), H(20), H(28)],
            vec![H(21), H(22), H(30)],
            vec![H(25), H(29), H(27), H(31)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
            vec![H(16), H(26), H(25)],
            vec![H(18), H(28), H(27)],
            vec![H(20), H(30), H(29)],
            vec![H(22), H(24), H(31)],
        ])).try_into_orientable();
        let (orientable, fans) = orientable.unwrap();
        assert_orientable_vh_isomorphic(&orientable, &ti_vec![
            vec![H(0), H(8) , H(7)],
            vec![H(2), H(10), H(1)],
            vec![H(4), H(12), H(3)],
            vec![H(6), H(14), H(5)],
            vec![H(9), H(11), H(13), H(15)],
            vec![H(16), H(24), H(23)],
            vec![H(18), H(26), H(17)],
            vec![H(20), H(28), H(19)],
            vec![H(22), H(30), H(21)],
            vec![H(25), H(27), H(29), H(31)],
        ], &ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
            vec![H(16), H(26), H(25)],
            vec![H(18), H(28), H(27)],
            vec![H(20), H(30), H(29)],
            vec![H(22), H(24), H(31)],
        ], &fans);

        // Three triangles from a single edge (not a orientable)
        let not_a_orientable = Frame {
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
        ])).try_into_orientable();
        assert!(matches!(not_a_orientable.unwrap_err()[0], NotAManifold(ManifoldError::EdgeAdjacentToTooManyFaceCorners { edge: E(0), face_corners: _ })));

        // A square cross with an edge sticking out from the middle. Not a orientable
        let not_a_orientable = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(13), H(11), H(15), H(16)],
            vec![H(17)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ])).try_into_orientable();
        assert!(matches!(not_a_orientable.unwrap_err()[0], NotAManifold(ManifoldError::VertexFacesDoNotFormFansOrLoop { vertex: V(4), face_corners: _ })));
    }
}