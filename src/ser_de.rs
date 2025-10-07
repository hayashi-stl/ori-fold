use std::{collections::HashMap, path::Path, sync::Arc};

use exact_number::{basis::ArcBasis, rat::{Integer, Natural, Rat}, BasedExpr};
use exact_number::malachite::base::num::basic::traits::{Zero, One};
use serde::{Deserialize, Serialize};
use serde::de::Error;
use serde_json::Value;

use crate::{fold::{EdgeAssignment, EdgeOrder, FaceOrder, FileClass, Fold, Frame, FrameAttribute, FrameClass, FrameUnit}, ser_de::basis::deserialize};

#[derive(Clone, Debug)]
struct SerDeRat(Rat);

impl Serialize for SerDeRat {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error> where S: serde::Serializer {
        let (numer, denom) = self.0.clone().into_signed_numerator_and_denominator();
        [numer.to_string(), denom.to_string()].serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for SerDeRat {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error> where D: serde::Deserializer<'de> {
        let fraction = String::deserialize(deserializer)?;
        let split = fraction.split("/").collect::<Vec<_>>();
        let numer = split[0].parse::<Integer>().map_err(|e| D::Error::custom(format!("can't parse numerator: {e:?}")))?;
        let denom = split.get(1)
            .map(|s| s.parse::<Natural>().map_err(|e| D::Error::custom(format!("can't parse denominator: {e:?}"))))
            .transpose()?
            .unwrap_or(Natural::ONE);
        if denom == Natural::ZERO {
            Err(D::Error::custom("denominator is 0"))?;
        }
        Ok(SerDeRat(Rat::from_integers(numer, denom.into())))
    }
}

fn version() -> String { "1.2".to_owned() }

#[derive(Clone, Debug)]
#[derive(Serialize, Deserialize)]
pub struct SerDeFold {
    /// The version of the FOLD spec that the file assumes.
    /// **Strongly recommended**, in case we ever have to make
    /// backward-incompatible changes.
    #[serde(default = "version")]
    #[serde(with = "crate::ser_de::file_spec")]
    pub file_spec: String,
    /// The software that created the file.
    /// **Recommended** for files output by computer software;
    /// less important for files made by hand.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_creator: Option<String>,
    /// The human author
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_author: Option<String>,
    /// A title for the entire file.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_title: Option<String>,
    /// A description of the entire file.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_description: Option<String>,
    /// A subjective interpretation about what the entire file represents.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub file_classes: Vec<FileClass>,
    /// The frames in the file
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub file_frames: Vec<Frame>,
}

#[derive(Clone, Debug)]
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct SerDeFrame {
    /// The human author
    #[serde(default, skip_serializing_if = "Option::is_none")]
    frame_author: Option<String>,
    /// A title for the frame.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    frame_title: Option<String>,
    /// A description of the frame.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    frame_description: Option<String>,
    /// A subjective interpretation about what the frame represents.
    #[serde(default)]
    frame_classes: Vec<FrameClass>,
    /// Attributes that objectively describe properties of the
    /// folded structure being represented.
    #[serde(default)]
    frame_attributes: Vec<FrameAttribute>,
    /// Physical or logical unit that all coordinates are relative to
    #[serde(default)]
    frame_unit: FrameUnit,

    /// The basis that exact coordinates uses.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(with = "crate::ser_de::basis")]
    #[serde(rename = "frame_exact:basis")]
    basis: Option<ArcBasis>,

    /// For each vertex, an array of approximate coordinates,
    /// such as `[x, y, z]` or `[x, y]` (where `z` is implicitly zero).
    /// In higher dimensions, all trailing unspecified coordinates are implicitly
    /// zero.  **Recommended** except for frames with attribute `Abstract`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "vertices_coords")]
    vertices_coords_f64: Option<Vec<Vec<f64>>>,
    /// For each vertex, an array of exact coordinates,
    /// such as `[x, y, z]` or `[x, y]` (where `z` is implicitly zero).
    /// In higher dimensions, all trailing unspecified coordinates are implicitly
    /// zero.  **Recommended** except for frames with attribute `Abstract`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "vertices_exact:coords")]
    vertices_coords_exact: Option<Vec<Vec<Vec<SerDeRat>>>>,
    /// For each vertex, an array of vertices (vertex IDs)
    /// that are adjacent along edges.  If the frame represents an orientable
    /// manifold or planar linkage, this list should be ordered counterclockwise
    /// around the vertex (possibly repeating a vertex more than once).
    /// If the frame is a nonorientable manifold, this list should be cyclically
    /// ordered around the vertex (possibly repeating a vertex).
    /// Otherwise, the order is arbitrary.
    /// **Recommended** in any frame lacking `edges_vertices` property
    /// (otherwise `vertices_vertices` can easily be computed from
    /// `edges_vertices` as needed).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    vertices_vertices: Option<Vec<Vec<usize>>>,
    /// For each vertex, an array of edge IDs for the edges
    /// incident to the vertex.  If the frame represents an orientable manifold,
    /// this list should be ordered counterclockwise around the vertex.
    /// If the frame is a nonorientable manifold, this list should be cyclically
    /// ordered around the vertex.
    /// In all cases, the linear order should match `vertices_vertices` if both
    /// are specified: `vertices_edges[v][i]` should be an edge connecting vertices
    /// `v` and `vertices_vertices[v][i]`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    vertices_edges: Option<Vec<Vec<usize>>>,
    /// For each vertex, an array of face IDs for the faces
    /// incident to the vertex, possibly including `None`s.
    /// If the frame represents a manifold, `vertices_faces` should align with
    /// `vertices_vertices` and/or `vertices_edges`:
    /// `vertices_faces[v][i]` should be either
    ///
    /// * the face containing vertices
    ///   `vertices_vertices[v][i]` and `vertices_vertices[v][(i+1)%d]` and
    ///   containing edges `vertices_edges[v][i]` and `vertices_edges[v][(i+1)%d]`,
    ///   where `d` is the degree of vertex `v`; or
    /// * `None` if such a face doesn't exist.
    ///
    /// If the frame represents an orientable manifold,
    /// this list should be ordered counterclockwise around the vertex
    /// (possibly repeating a face more than once).  If the frame is a
    /// nonorientable manifold, this list should be cyclically ordered around the
    /// vertex (possibly repeating a vertex), and matching the cyclic order of
    /// `vertices_vertices` and/or `vertices_edges` (if either is specified).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    vertices_faces: Option<Vec<Vec<Option<usize>>>>,
    /// `edges_vertices`: For each edge, an array `[u, v]` of two vertex IDs for
    /// the two endpoints of the edge.  This effectively defines the *orientation*
    /// of the edge, from `u` to `v`.  (This orientation choice is arbitrary,
    /// but is used to define the ordering of `edges_faces`.)
    /// **Recommended** in frames having any `edges_...` property
    /// (e.g., to represent mountain-valley assignment).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    edges_vertices: Option<Vec<[usize; 2]>>,
    /// For each edge, an array of face IDs for the faces incident
    /// to the edge, possibly including `None`s.
    /// For nonmanifolds in particular, the `Some()` faces should be listed in
    /// counterclockwise order around the edge,
    /// relative to the orientation of the edge.
    /// For manifolds, the array for each edge should be an array of length 2,
    /// where the first entry is the face locally to the "left" of the edge
    /// (or `None` if there is no such face) and the second entry is the face
    /// locally to the "right" of the edge (or `None` if there is no such face);
    /// for orientable manifolds, "left" and "right" must be consistent with the
    /// manifold orientation given by the counterclockwise orientation of faces.
    /// However, a boundary edge may also be represented by a length-1 array, with
    /// the `None` omitted, to be consistent with the nonmanifold representation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    edges_faces: Option<Vec<Vec<Option<usize>>>>,
    /// For each edge, a string representing its fold direction assignment
    #[serde(default, skip_serializing_if = "Option::is_none")]
    edges_assignment: Option<Vec<EdgeAssignment>>,
    /// For each edge, the fold angle (deviation from flatness)
    /// along each edge of the pattern.  The fold angle is a number in degrees
    /// lying in the range [-180, 180].  The fold angle is positive for
    /// valley folds, negative for mountain folds, and zero for flat, unassigned,
    /// and border folds.  Accordingly, the sign of `edge_foldAngle` should match
    /// `edges_assignment` if both are specified.
    /// 
    /// For now, this is an `f64` regardless of `N`, because finding a good representation
    /// with `N` is tricky. `[cos(a), sin(a)]` doesn't work because -180° and 180° are both in range.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edges_foldAngle")]
    edges_fold_angle_f64: Option<Vec<f64>>,
    /// For each edge, the *exact* fold angle (deviation from flatness)
    /// along each edge of the pattern, written as `(negative?, cos(a), sin(a))`.
    /// Both `cos(a)` and `sin(a)` are stored to ensure they're both in `N`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edges_exact:foldAngle")]
    edges_fold_angle_cs: Option<Vec<(bool, Vec<SerDeRat>, Vec<SerDeRat>)>>,

    /// For each edge, the approximate length of the edge.
    /// This is mainly useful for defining the intrinsic geometry of
    /// abstract complexes where `vertices_coords` are unspecified;
    /// otherwise, `edges_length` can be computed from `vertices_coords`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edges_lengths")]
    edges_length_f64: Option<Vec<f64>>,
    /// For each edge, the exact squared length of the edge.
    /// This is mainly useful for defining the intrinsic geometry of
    /// abstract complexes where `vertices_coords` are unspecified;
    /// otherwise, `edges_length` can be computed from `vertices_coords`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edges_exact:lengths")]
    edges_length2_exact: Option<Vec<Vec<SerDeRat>>>,

    /// For each face, an array of vertex IDs for the vertices
    /// around the face *in counterclockwise order*.  This array can repeat the
    /// same vertex multiple times (e.g., if the face has a "slit" in it).
    /// **Recommended** in any frame having faces.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    faces_vertices: Option<Vec<Vec<usize>>>,
    /// For each face, an array of edge IDs for the edges around
    /// the face *in counterclockwise order*.  In addition to the matching cyclic
    /// order, `faces_vertices` and `faces_edges` should align in start so that
    /// `faces_edges[f][i]` is the edge connecting `faces_vertices[f][i]` and
    /// `faces_vertices[f][(i+1)%d]` where `d` is the degree of face `f`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    faces_edges: Option<Vec<Vec<usize>>>,
    /// `faces_faces`: For each face, an array of face IDs for the faces *sharing
    /// edges* around the face, possibly including `null`s.
    /// If the frame is a manifold, the faces should be listed in counterclockwise
    /// order and in the same linear order as `faces_edges` (if it is specified):
    /// `f` and `faces_faces[f][i]` should be the faces incident to the edge
    /// `faces_edges[f][i]`, unless that edge has no face on the other side,
    /// in which case `faces_faces[f][i]` should be `null`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    faces_faces: Option<Vec<Vec<Option<usize>>>>,
    
    /// An array of triples `(f, g, s)` where `f` and `g` are face IDs
    /// and `s` is a `FaceOrder`:
    /// * `Above` indicates that face `f` lies *above* face `g`,
    ///   i.e., on the side pointed to by `g`'s normal vector in the folded state.
    /// * `Below` indicates that face `f` lies *below* face `g`,
    ///   i.e., on the side opposite `g`'s normal vector in the folded state.
    /// * `Unknown` indicates that `f` and `g` have unknown stacking order
    ///   (e.g., they do not overlap in their interiors).
    ///
    /// **Recommended** for frames with interior-overlapping faces.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    face_orders: Option<Vec<(usize, usize, FaceOrder)>>,
    /// An array of triples `[e, f, s]` where `e` and `f` are edge IDs
    /// and `s` is a `EdgeOrder`.
    /// * `Left` indicates that edge `e` lies locally on the *left* side of edge `f`
    ///   (relative to edge `f`'s orientation given by `edges_vertices`)
    /// * `Right` indicates that edge `e` lies locally on the *right* side of edge
    ///   `f` (relative to edge `f`'s orientation given by `edges_vertices`)
    /// * 0 indicates that `e` and `f` have unknown stacking order
    ///   (e.g., they do not overlap in their interiors).
    ///
    /// This property makes sense only in 2D.
    /// **Recommended** for linkage configurations with interior-overlapping edges.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    edge_orders: Option<Vec<(usize, usize, EdgeOrder)>>,

    /// Parent frame ID.  Intuitively, this frame (the child)
    /// is a modification (or, in general, is related to) the parent frame.
    /// This property is optional, but enables organizing frames into a tree
    /// structure.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    frame_parent: Option<usize>,
    /// If true, any properties in the parent frame
    /// (or recursively inherited from an ancestor) that is not overridden in
    /// this frame are automatically inherited, allowing you to avoid duplicated
    /// data in many cases.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    frame_inherit: Option<bool>,
    /// Custom data
    #[serde(flatten)]
    frame_custom: HashMap<String, Value>,
}

fn based_expr_from_de(coeffs: Vec<SerDeRat>, basis: ArcBasis) -> Result<BasedExpr, String> {
    let coeffs = coeffs.into_iter().map(|s| s.0).collect::<Vec<_>>();
    BasedExpr::try_from_coeffs_and_basis(coeffs, basis).ok_or("too many coefficients for the basis".to_owned())
}

impl TryFrom<SerDeFrame> for Frame {
    type Error = String;

    fn try_from(value: SerDeFrame) -> Result<Self, Self::Error> {
        if (value.vertices_coords_exact.is_some() || value.edges_fold_angle_cs.is_some() || value.edges_length2_exact.is_some())
            && value.basis.is_none()
        {
            Err("exact coordinates without a basis")?
        }
        Ok(Self {
            frame_author: value.frame_author,
            frame_title: value.frame_title,
            frame_description: value.frame_description,
            frame_classes: value.frame_classes,
            frame_attributes: value.frame_attributes,
            frame_unit: value.frame_unit,
            basis: value.basis.clone(),
            vertices_coords_f64: value.vertices_coords_f64,
            vertices_coords_exact: value.vertices_coords_exact
                .map(|vs| vs.into_iter()
                    .map(|v| v.into_iter()
                        .map(|c| based_expr_from_de(c, value.basis.as_ref().unwrap().clone()))
                        .collect::<Result<Vec<_>, _>>())
                    .collect::<Result<Vec<_>, _>>())
                .transpose()?,
            vertices_vertices: value.vertices_vertices,
            vertices_edges: value.vertices_edges,
            vertices_faces: value.vertices_faces,
            edges_vertices: value.edges_vertices,
            edges_faces: value.edges_faces,
            edges_assignment: value.edges_assignment,
            edges_fold_angle_f64: value.edges_fold_angle_f64,
            edges_length_f64: value.edges_length_f64,
            edges_fold_angle_cs: value.edges_fold_angle_cs
                .map(|ls| ls.into_iter()
                    .map(|(sign, cos, sin)|
                        based_expr_from_de(cos, value.basis.as_ref().unwrap().clone())
                            .and_then(|cos| based_expr_from_de(sin, value.basis.as_ref().unwrap().clone()).map(|sin| (cos, sin)))
                            .map(|(cos, sin)| (sign, cos, sin)))
                    .collect::<Result<Vec<_>, _>>())
                .transpose()?,
            edges_length2_exact: value.edges_length2_exact
                .map(|ls| ls.into_iter()
                    .map(|c| based_expr_from_de(c, value.basis.as_ref().unwrap().clone()))
                    .collect::<Result<Vec<_>, _>>())
                .transpose()?,
            faces_vertices: value.faces_vertices,
            faces_edges: value.faces_edges,
            faces_faces: value.faces_faces,
            face_orders: value.face_orders,
            edge_orders: value.edge_orders,
            frame_parent: value.frame_parent,
            frame_inherit: value.frame_inherit,
            frame_custom: value.frame_custom,
        })
    }
}

fn based_expr_to_ser(b: BasedExpr) -> Vec<SerDeRat> {
    match b {
        BasedExpr::Baseless(q) => vec![SerDeRat(q)],
        // into_iter doesn't do the expected thing
        BasedExpr::Based(mut coeffs, _) => coeffs.iter_mut().map(|c| std::mem::take(c)).map(SerDeRat).collect::<Vec<_>>(),
    }
}

impl From<Frame> for SerDeFrame {
    fn from(value: Frame) -> Self {
        Self {
            frame_author: value.frame_author,
            frame_title: value.frame_title,
            frame_description: value.frame_description,
            frame_classes: value.frame_classes,
            frame_attributes: value.frame_attributes,
            frame_unit: value.frame_unit,
            basis: value.basis,
            vertices_coords_f64: value.vertices_coords_f64,
            vertices_coords_exact: value.vertices_coords_exact
                .map(|vs| vs.into_iter()
                    .map(|v| v.into_iter().map(based_expr_to_ser).collect::<Vec<_>>())
                    .collect::<Vec<_>>()),
            vertices_vertices: value.vertices_vertices,
            vertices_edges: value.vertices_edges,
            vertices_faces: value.vertices_faces,
            edges_vertices: value.edges_vertices,
            edges_faces: value.edges_faces,
            edges_assignment: value.edges_assignment,
            edges_fold_angle_f64: value.edges_fold_angle_f64,
            edges_length_f64: value.edges_length_f64,
            edges_fold_angle_cs: value.edges_fold_angle_cs
                .map(|ls| ls.into_iter()
                    .map(|(sign, cos, sin)|
                        (sign, based_expr_to_ser(cos), based_expr_to_ser(sin)))
                    .collect::<Vec<_>>()),
            edges_length2_exact: value.edges_length2_exact
                .map(|ls| ls.into_iter().map(based_expr_to_ser).collect::<Vec<_>>()),
            faces_vertices: value.faces_vertices,
            faces_edges: value.faces_edges,
            faces_faces: value.faces_faces,
            face_orders: value.face_orders,
            edge_orders: value.edge_orders,
            frame_parent: value.frame_parent,
            frame_inherit: value.frame_inherit,
            frame_custom: value.frame_custom,
        }
    }
}

impl Fold {
    /// Reads this FOLD from a file.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, serde_json::error::Error> {
        use serde_json::error::Error;
        use serde::de::{Error as _};

        let file = std::fs::read_to_string(path.as_ref())
            .map_err(|e| Error::custom(format!("cannot read from {:?}: {:?}", path.as_ref(), e)))?;
        Fold::from_str(&file)
    }

    /// #[serde(flatten)] is glitched for structs, so we have to work around it
    pub fn from_str(s: &str) -> Result<Self, serde_json::error::Error> {
        // with a double deserialize!
        let fold = serde_json::de::from_str::<SerDeFold>(s)?;
        let key_frame = serde_json::de::from_str::<Frame>(s)?;
        let mut fold = Fold {
            file_spec: fold.file_spec,
            file_creator: fold.file_creator,
            file_author: fold.file_author,
            file_title: fold.file_title,
            file_description: fold.file_description,
            file_classes: fold.file_classes,
            key_frame: key_frame,
            file_frames: fold.file_frames,
            file_custom: HashMap::new()
        };

        // Move custom fields appropriately
        let mut file_custom = HashMap::new();
        let mut frame_custom = HashMap::new();
        for keys in std::mem::take(&mut fold.key_frame.frame_custom) {
            if keys.0.starts_with("frame_") {
                frame_custom.insert(keys.0, keys.1);
            } else if keys.0.contains(":") {
                file_custom.insert(keys.0, keys.1);
            }
        }
        
        fold.file_custom = file_custom;
        fold.key_frame.frame_custom = frame_custom;
        Ok(fold)
    }
}

/// For whatever reason, the `file_spec` is specified to be a number and not a string, requiring custom (de)serialization
pub mod file_spec {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use serde::ser::{Error as _};

    pub fn serialize<S>(spec: &String, ser: S) -> Result<S::Ok, S::Error> where S: Serializer {
        let spec = spec.parse::<f64>()
            .map_err(|e| S::Error::custom(format!("unfortunately spec must be convertible to a number: {e:?}")))?;
        spec.serialize(ser)
    }

    pub fn deserialize<'de, D>(de: D) -> Result<String, D::Error> where D: Deserializer<'de> {
        println!("Deserializing!");
        Ok(f64::deserialize(de)?.to_string())
    }
}

pub mod basis {
    use exact_number::basis::{ArcBasis, Basis, SqrtExpr};
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use serde::de::Error;

    pub fn serialize<S>(basis: &Option<ArcBasis>, ser: S) -> Result<S::Ok, S::Error> where S: Serializer {
        let basis = basis.as_ref();
        let basis = if let Some(b) = basis { b } else {
            return ser.serialize_none();
        };
        ser.serialize_some(
            &basis.exprs.iter().map(|expr| expr.to_string()).collect::<Vec<_>>()
        )
    }

    pub fn deserialize<'de, D>(de: D) -> Result<Option<ArcBasis>, D::Error> where D: Deserializer<'de> {
        let strings = Vec::<String>::deserialize(de)?;
        let exprs = strings.into_iter().map(|b| b.parse::<SqrtExpr>())
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| D::Error::custom(format!("failed to parse SqrtExpr: {e:?}")))?;
        Basis::new_arc_checked(exprs)
            .map(Some)
            .map_err(|e| D::Error::custom(format!("not a basis: {e:?}")))
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;

    use exact_number::{based_expr, basis::{ArcBasis, Basis}, sqrt_expr};
    use serde_json::Value;

    use crate::fold::{FileClass, Fold};

    #[test]
    fn test_load_partial_metadata() {
        let fold = Fold::from_file("test/x.fold").unwrap();
        assert_eq!(fold.file_spec, "1.1");
        assert_eq!(fold.file_creator, Some("oriedita".to_owned()));
        assert_eq!(fold.file_author, None);
        assert_eq!(fold.file_title, None);
        assert_eq!(fold.file_description, None);
        assert_eq!(fold.file_classes, vec![]);
        assert_eq!(fold.file_custom, vec![("oriedita:version".to_owned(), Value::String("1.1.3".to_owned()))].into_iter().collect::<HashMap<_, _>>());
    }

    #[test]
    fn test_load_full_metadata() {
        let fold = Fold::from_file("test/x-all-file-metadata.fold").unwrap();
        assert_eq!(fold.file_spec, "1.1");
        assert_eq!(fold.file_creator, Some("oriedita".to_owned()));
        assert_eq!(fold.file_author, Some("hayastl".to_owned()));
        assert_eq!(fold.file_title, Some("x".to_owned()));
        assert_eq!(fold.file_description, Some("a simple cross pattern".to_owned()));
        assert_eq!(fold.file_classes, vec![FileClass::SingleModel]);
        assert_eq!(fold.file_custom, vec![("oriedita:version".to_owned(), Value::String("1.1.3".to_owned()))].into_iter().collect::<HashMap<_, _>>());
    }

    #[test]
    fn test_load_exact_coords() {
        let fold = Fold::from_file("test/x-exact-coords.fold").unwrap();
        let basis = Basis::new_arc(vec![sqrt_expr!(1)]);
        let coords = vec![
            vec![based_expr!(-200), based_expr!(-200)],
            vec![based_expr!(-200), based_expr!(200)],
            vec![based_expr!(200), based_expr!(-200)],
            vec![based_expr!(200), based_expr!(200)],
            vec![based_expr!(0), based_expr!(0)],
        ];
        assert_eq!(fold.key_frame.basis, Some(basis));
        assert_eq!(fold.key_frame.vertices_coords_exact, Some(coords));
    }

    #[test]
    fn test_load_exact_coords_bird_base() {
        let fold = Fold::from_file("test/x-exact-coords-bird-base.fold").unwrap();
        let basis = Basis::new_arc(vec![sqrt_expr!(1), sqrt_expr!(2)]);
        let coords = vec![
            vec![based_expr!(0 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2)],
            vec![based_expr!(1 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2)],
            vec![based_expr!(1 + 0 sqrt 2), based_expr!(1 + 0 sqrt 2)],
            vec![based_expr!(0 + 0 sqrt 2), based_expr!(1 + 0 sqrt 2)],
            vec![based_expr!(1/2 + 0 sqrt 2), based_expr!(1/2 + 0 sqrt 2)],
            vec![based_expr!(1/2 + 0 sqrt 2), based_expr!(-1/2 + 1/2 sqrt 2)],
            vec![based_expr!(3/2 - 1/2 sqrt 2), based_expr!(1/2 + 0 sqrt 2)],
            vec![based_expr!(1/2 + 0 sqrt 2), based_expr!(3/2 - 1/2 sqrt 2)],
            vec![based_expr!(-1/2 + 1/2 sqrt 2), based_expr!(1/2 + 0 sqrt 2)],
        ];
        assert_eq!(fold.key_frame.basis, Some(basis));
        assert_eq!(fold.key_frame.vertices_coords_exact, Some(coords));
    }

    #[test]
    fn test_load_invalid_exact_coords_baseless() {
        let fold = Fold::from_file("test/invalid-x-exact-coords-baseless.fold");
        assert!(fold.is_err());
    }

    #[test]
    fn test_load_invalid_empty() {
        let fold = Fold::from_file("test/invalid-x-empty.fold");
        assert!(fold.is_err());
    }
}