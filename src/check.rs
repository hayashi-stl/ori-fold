//! For checking FOLD data against the spec.
//! There's a bunch of conditions, especially for manifold FOLDs

use nalgebra::{DMatrix, DMatrixView, DVector};

use crate::fold::{Fold, Frame, FrameAttribute};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum CheckError {
    FrameAttributeConflict{ conflict: [FrameAttribute; 2] },
    MissingCoordinates{ attribute: FrameAttribute },
    DimensionMismatch{ attribute: FrameAttribute, dims: usize }
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

    /// Check the frame attributes and get a list of the errors.
    /// 
    /// If this function returns no errors, you can assume that
    /// if `FrameAttribute::_2D` or `FrameAttribute::_3D` are specified,
    /// then there are vertex coordinates with the corresponding dimension.
    pub fn frame_attribute_vertex_coords_conflicts(&self) -> Vec<CheckError> {
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

    /// Check the frame and get a list of the errors.
    pub fn check(&self) -> Vec<CheckError> {
        let mut errors = vec![];
        errors.extend(self.frame_attribute_conflicts());
        errors.extend(self.frame_attribute_vertex_coords_conflicts());
        errors
    }
}

#[cfg(test)]
mod test {
    use exact_number::based_expr;
    use nalgebra::DMatrix;

    use crate::{check::CheckError, fold::{Frame, FrameAttribute}};

    fn frame_with_attributes(attributes: impl IntoIterator<Item = FrameAttribute>) -> Frame {
        Frame {
            frame_attributes: attributes.into_iter().collect::<Vec<_>>(),
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
    fn test_check_frame_attribute_vertex_coords_conflicts() {
        use FrameAttribute::*;
        use CheckError::*;

        assert_eq!(canonicalize(frame_with_attributes([_2D]).frame_attribute_vertex_coords_conflicts()),
            canonicalize(vec![MissingCoordinates { attribute: FrameAttribute::_2D }]));
        assert_eq!(canonicalize(frame_with_attributes([_3D]).frame_attribute_vertex_coords_conflicts()),
            canonicalize(vec![MissingCoordinates { attribute: FrameAttribute::_3D }]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: vec![_2D],
                vertices_coords_exact: Some(DMatrix::from_vec(3, 1, vec![based_expr!(1), based_expr!(2), based_expr!(3)])),
                ..Default::default()
            }.frame_attribute_vertex_coords_conflicts()),
            canonicalize(vec![DimensionMismatch { attribute: _2D, dims: 3 }]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: vec![_3D],
                vertices_coords_exact: Some(DMatrix::from_vec(2, 2, vec![based_expr!(1), based_expr!(2), based_expr!(3), based_expr!(4)])),
                ..Default::default()
            }.frame_attribute_vertex_coords_conflicts()),
            canonicalize(vec![DimensionMismatch { attribute: _3D, dims: 2 }]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: vec![_3D],
                vertices_coords_exact: Some(DMatrix::from_vec(3, 1, vec![based_expr!(1), based_expr!(2), based_expr!(3)])),
                ..Default::default()
            }.frame_attribute_vertex_coords_conflicts()),
            canonicalize(vec![]));
        assert_eq!(canonicalize(Frame {
                frame_attributes: vec![_2D],
                vertices_coords_exact: Some(DMatrix::from_vec(2, 2, vec![based_expr!(1), based_expr!(2), based_expr!(3), based_expr!(4)])),
                ..Default::default()
            }.frame_attribute_vertex_coords_conflicts()),
            canonicalize(vec![]));
    }
}