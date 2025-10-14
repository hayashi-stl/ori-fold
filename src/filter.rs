use crate::fold::{CoordsRef, Frame};

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