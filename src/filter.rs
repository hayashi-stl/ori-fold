use crate::fold::{CoordsRef, Edge, Frame, HalfEdge, Vertex};

impl Frame {
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
    use exact_number::based_expr;
    use indexmap::IndexMap;
    use nalgebra::DMatrix;
    use serde_json::json;

    use crate::{filter::{edges_vertices_incident, other_vertex}, fold::Frame};
    use crate::fold::{Vertex as V, Edge as E, Face as F};

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
}