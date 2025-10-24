//! For checking FOLD data against the spec.
//! There's a bunch of conditions, especially for manifold FOLDs

use std::fmt::Display;

use exact_number::{BasedExpr, Angle};
use nalgebra::{DMatrix, DMatrixView, DVector};
use typed_index_collections::TiSlice;

use crate::{fold::{Edge, EdgeAssignment, Fold, Frame, FrameAttribute}, manifold::{ManifoldError, OrientableError}};

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(PartialEq, Eq, PartialOrd, Ord))] // Really no point in this derivation outside of tests
pub enum CheckError {
    FrameAttributeConflict{ conflict: [FrameAttribute; 2] },
    MissingCoordinates{ attribute: FrameAttribute },
    DimensionMismatch{ attribute: FrameAttribute, dims: usize },
    NotAManifold(ManifoldError),
    NotOrientable(OrientableError),
    InvalidEdgeAssignment{ attribute: FrameAttribute, edge: Edge, assignment: EdgeAssignment },
    ExactFoldAngleIsNotNormalized{ edge: Edge, angle: Angle },
}

impl Display for CheckError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::FrameAttributeConflict { conflict } =>
                write!(f, "the following attributes conflict: {conflict:?}"),
            Self::MissingCoordinates { attribute } =>
                write!(f, "the attribute {attribute:?} is specified but vertex coordinates are missing"),
            Self::DimensionMismatch { attribute, dims } =>
                write!(f, "the attribute {attribute:?} is specified but vertex coordinates are {dims}D"),
            Self::NotAManifold(error) => write!(f, "{error}"),
            Self::NotOrientable(error) => write!(f, "{error}"),
            Self::InvalidEdgeAssignment { attribute, edge, assignment } =>
                write!(f, "the attribute {attribute:?} is specified but edge {edge} has assignment {assignment:?}"),
            Self::ExactFoldAngleIsNotNormalized { edge, angle } =>
                write!(f, "the exact angle of {edge} is not normalized: {angle:?}"),
        }
    }
}

impl Frame {
    /// Check the frame attributes and get a list of the errors.
    /// 
    /// If this function returns no errors, you can assume that there
    /// are no conflicting attributes in this frame.
    pub fn frame_attribute_conflicts(&self) -> Vec<CheckError> {
        let conflicts = [
            [FrameAttribute::_2D, FrameAttribute::_3D],
            [FrameAttribute::_2D, FrameAttribute::Abstract],
            [FrameAttribute::_3D, FrameAttribute::Abstract],
            [FrameAttribute::Manifold, FrameAttribute::NonManifold],
            [FrameAttribute::Orientable, FrameAttribute::NonManifold],
            [FrameAttribute::Orientable, FrameAttribute::NonOrientable],
            [FrameAttribute::SelfTouching, FrameAttribute::NonSelfTouching],
            [FrameAttribute::SelfIntersecting, FrameAttribute::NonSelfTouching],
            [FrameAttribute::SelfIntersecting, FrameAttribute::NonSelfIntersecting],
            [FrameAttribute::Cuts, FrameAttribute::NoCuts],
            [FrameAttribute::Joins, FrameAttribute::NoJoins],
            [FrameAttribute::ConvexFaces, FrameAttribute::NonConvexFaces],
        ];
    
        let mut errors = vec![];
        for conflict in &conflicts {
            if self.frame_attributes.contains(&conflict[0]) && self.frame_attributes.contains(&conflict[1]) {
                errors.push(CheckError::FrameAttributeConflict { conflict: conflict.clone() });
            }
        }
        errors
    }

    /// Check the dimensional frame attributes and get a list of the errors.
    /// 
    /// If this function returns no errors, you can assume that
    /// if `FrameAttribute::_2D` or `FrameAttribute::_3D` are specified,
    /// then there are vertex coordinates with the corresponding dimension.
    pub fn frame_attribute_dimensions_conflicts(&self) -> Vec<CheckError> {
        let vertex_coords = if let Some(c) = self.coords_ref() { c } else {
            return if self.frame_attributes.contains(&FrameAttribute::_2D) {
                vec![CheckError::MissingCoordinates { attribute: FrameAttribute::_2D }]
            } else if self.frame_attributes.contains(&FrameAttribute::_3D) {
                vec![CheckError::MissingCoordinates { attribute: FrameAttribute::_3D }]
            } else {
                vec![]
            }
        };
        let dims = [
            (FrameAttribute::_2D, 2),
            (FrameAttribute::_3D, 3)
        ];
        let mut errors = vec![];
        for (attribute, dims) in dims {
            if self.frame_attributes.contains(&attribute) && vertex_coords.num_dimensions() != dims {
                errors.push(CheckError::DimensionMismatch { attribute, dims: vertex_coords.num_dimensions() })
            }
        }
        errors
    }

    /// Check the cuts/joins frame attributes and get a list of the errors.
    /// 
    /// If this function returns no errors, you can assume that
    /// if `FrameAttribute::NoCuts` or `FrameAttribute::NoJoins` is specified,
    /// then there are in fact no cut edges or no join edges, respectively.
    pub fn frame_attribute_edge_assignment_conflicts(&self) -> Vec<CheckError> {
        let conflicts = [
            (FrameAttribute::NoCuts, EdgeAssignment::Cut),
            (FrameAttribute::NoJoins, EdgeAssignment::Join),
        ];

        let edges_assignment = if let Some(ea) = self.edges_assignment.as_ref() { ea } else { return vec![] };
        let mut errors = vec![];
        for (attr, bad_assigment) in conflicts {
            if !self.frame_attributes.contains(&attr) { continue }
            
            for (e, &assignment) in edges_assignment.iter_enumerated() {
                if assignment == bad_assigment {
                    errors.push(CheckError::InvalidEdgeAssignment { attribute: attr.clone(), edge: e, assignment: assignment });
                }
            }
        }

        errors
    }

    /// Check the frame and get a list of the errors.
    /// Modifies the representation as necessary to handle frame attributes.
    pub fn check(mut self) -> Result<Self, Vec<CheckError>> {
        let mut errors = vec![];
        errors.extend(self.frame_attribute_conflicts());
        errors.extend(self.frame_attribute_dimensions_conflicts());
        errors.extend(self.frame_attribute_edge_assignment_conflicts());

        if self.remove_attribute(FrameAttribute::Orientable) {
            self = self.try_into_orientable()
                .map(|m| m.0)
                .map_err(|e| errors.clone().into_iter().chain(e.into_iter().map(CheckError::NotOrientable)).collect::<Vec<_>>())?;
        } else if self.remove_attribute(FrameAttribute::Manifold) {
            self = self.try_into_manifold()
                .map(|m| m.0)
                .map_err(|e| errors.clone().into_iter().chain(e.into_iter().map(CheckError::NotAManifold)).collect::<Vec<_>>())?;
        }

        // All other frame attribute checks require geometry

        // TODO: Check for self-intersection. This is complicated, especially when self-touching is allowed.
        // TODO: Check face_orders and edge_orders for consistency. This is a little complicated too.
        
        if errors.is_empty() { Ok(self) } else { Err(errors) }
    }
}

#[cfg(test)]
mod test {
    use exact_number::based_expr;
    use indexmap::{indexset, IndexSet};
    use nalgebra::DMatrix;
    use typed_index_collections::ti_vec;

    use crate::{check::CheckError, fold::{EdgeAssignment, Frame, FrameAttribute}};
    use crate::fold::{Edge as E};

    fn frame_with_attributes(attributes: impl IntoIterator<Item = FrameAttribute>) -> Frame {
        Frame {
            frame_attributes: attributes.into_iter().collect::<IndexSet<_>>(),
            ..Default::default()
        }
    }

    fn canonicalize(mut errors: Vec<CheckError>) -> Vec<CheckError> {
        use CheckError::*;

        for error in &mut errors {
            match error {
                FrameAttributeConflict { conflict } => conflict.sort(),
                MissingCoordinates { .. } => {},
                DimensionMismatch { .. } => {},
                _ => {},
            }
        }
        errors.sort();
        errors
    }

    #[test]
    fn test_check_frame_attribute_conflicts() {
        use FrameAttribute::*;
        use CheckError::*;

        assert_eq!(canonicalize(frame_with_attributes([_2D, Orientable, Manifold, NonConvexFaces]).frame_attribute_conflicts()),
            canonicalize(vec![]));
        assert_eq!(canonicalize(frame_with_attributes([_2D, Abstract]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [_2D, Abstract] }]));
        assert_eq!(canonicalize(frame_with_attributes([_3D, Abstract]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [_3D, Abstract] }]));
        assert_eq!(canonicalize(frame_with_attributes([_3D, _2D]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [_3D, _2D] }]));
        assert_eq!(canonicalize(frame_with_attributes([Manifold, NonManifold]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [Manifold, NonManifold] }]));
        assert_eq!(canonicalize(frame_with_attributes([Orientable, NonManifold]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [Orientable, NonManifold] }]));
        assert_eq!(canonicalize(frame_with_attributes([Orientable, NonOrientable]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [Orientable, NonOrientable] }]));
        assert_eq!(canonicalize(frame_with_attributes([SelfTouching, NonSelfTouching]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [SelfTouching, NonSelfTouching] }]));
        assert_eq!(canonicalize(frame_with_attributes([SelfIntersecting, NonSelfTouching]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [SelfIntersecting, NonSelfTouching] }]));
        assert_eq!(canonicalize(frame_with_attributes([SelfIntersecting, NonSelfIntersecting]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [SelfIntersecting, NonSelfIntersecting] }]));
        assert_eq!(canonicalize(frame_with_attributes([Cuts, NoCuts]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [Cuts, NoCuts] }]));
        assert_eq!(canonicalize(frame_with_attributes([Joins, NoJoins]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [Joins, NoJoins] }]));
        assert_eq!(canonicalize(frame_with_attributes([NonConvexFaces, ConvexFaces]).frame_attribute_conflicts()),
            canonicalize(vec![FrameAttributeConflict { conflict: [NonConvexFaces, ConvexFaces] }]));
        assert_eq!(canonicalize(frame_with_attributes([_2D, _3D, Abstract]).frame_attribute_conflicts()),
            canonicalize(vec![
                FrameAttributeConflict { conflict: [_2D, Abstract] },
                FrameAttributeConflict { conflict: [_3D, Abstract] },
                FrameAttributeConflict { conflict: [_3D, _2D] },
                ]));
    }

    #[test]
    fn test_check_frame_attribute_dimensions_conflicts() {
        use FrameAttribute::*;
        use CheckError::*;

        assert_eq!(canonicalize(frame_with_attributes([_2D]).frame_attribute_dimensions_conflicts()),
            canonicalize(vec![MissingCoordinates { attribute: FrameAttribute::_2D }]));
        assert_eq!(canonicalize(frame_with_attributes([_3D]).frame_attribute_dimensions_conflicts()),
            canonicalize(vec![MissingCoordinates { attribute: FrameAttribute::_3D }]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![_2D],
                vertices_coords_exact: Some(DMatrix::from_vec(3, 1, vec![based_expr!(1), based_expr!(2), based_expr!(3)])),
                ..Default::default()
            }.frame_attribute_dimensions_conflicts()),
            canonicalize(vec![DimensionMismatch { attribute: _2D, dims: 3 }]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![_3D],
                vertices_coords_exact: Some(DMatrix::from_vec(2, 2, vec![based_expr!(1), based_expr!(2), based_expr!(3), based_expr!(4)])),
                ..Default::default()
            }.frame_attribute_dimensions_conflicts()),
            canonicalize(vec![DimensionMismatch { attribute: _3D, dims: 2 }]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![_3D],
                vertices_coords_exact: Some(DMatrix::from_vec(3, 1, vec![based_expr!(1), based_expr!(2), based_expr!(3)])),
                ..Default::default()
            }.frame_attribute_dimensions_conflicts()),
            canonicalize(vec![]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![_2D],
                vertices_coords_exact: Some(DMatrix::from_vec(2, 2, vec![based_expr!(1), based_expr!(2), based_expr!(3), based_expr!(4)])),
                ..Default::default()
            }.frame_attribute_dimensions_conflicts()),
            canonicalize(vec![]));
    }

    #[test]
    fn test_check_frame_attribute_edge_assignment_conflicts() {
        use FrameAttribute::*;
        use EdgeAssignment::*;
        use CheckError::*;

        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoCuts],
                edges_assignment: Some(ti_vec![
                    Boundary
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoCuts],
                edges_assignment: Some(ti_vec![
                    Cut
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![
                InvalidEdgeAssignment { attribute: NoCuts, edge: E(0), assignment: Cut }
            ]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoCuts],
                edges_assignment: Some(ti_vec![
                    Cut,
                    Boundary,
                    Cut
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![
                InvalidEdgeAssignment { attribute: NoCuts, edge: E(0), assignment: Cut },
                InvalidEdgeAssignment { attribute: NoCuts, edge: E(2), assignment: Cut },
            ]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoJoins],
                edges_assignment: Some(ti_vec![
                    Boundary
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoJoins],
                edges_assignment: Some(ti_vec![
                    Join
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![
                InvalidEdgeAssignment { attribute: NoJoins, edge: E(0), assignment: Join }
            ]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoJoins],
                edges_assignment: Some(ti_vec![
                    Join,
                    Boundary,
                    Join
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![
                InvalidEdgeAssignment { attribute: NoJoins, edge: E(0), assignment: Join },
                InvalidEdgeAssignment { attribute: NoJoins, edge: E(2), assignment: Join },
            ]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoCuts],
                edges_assignment: Some(ti_vec![
                    Join,
                    Boundary,
                    Join
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoJoins],
                edges_assignment: Some(ti_vec![
                    Cut,
                    Boundary,
                    Cut
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: indexset![NoJoins, NoCuts],
                edges_assignment: Some(ti_vec![
                    Join,
                    Boundary,
                    Join,
                    Cut,
                ]),
                ..Default::default()
            }.frame_attribute_edge_assignment_conflicts()),
            canonicalize(vec![
                InvalidEdgeAssignment { attribute: NoCuts,  edge: E(3), assignment: Cut  },
                InvalidEdgeAssignment { attribute: NoJoins, edge: E(0), assignment: Join },
                InvalidEdgeAssignment { attribute: NoJoins, edge: E(2), assignment: Join },
            ]));
    }
}