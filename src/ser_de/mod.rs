use std::{path::Path, sync::Arc};

use exact_number::{basis::ArcBasis, rat::{Integer, Natural, Rat}, BasedExpr, Angle};
use exact_number::malachite::base::num::basic::traits::{Zero, One};
use indexmap::{IndexMap, IndexSet};
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use serde::de::Error;
use serde_json::{Map, Value};
use typed_index_collections::TiVec;

use crate::{fold::{Edge, EdgeAssignment, EdgeOrder, Face, FaceOrder, FileClass, Fold, Frame, FrameAttribute, FrameClass, FrameIndex, FrameUnit, Vertex}, ser_de::basis::deserialize};

mod pretty;
mod convert;

#[derive(Clone, Debug)]
pub(crate) struct SerDeRat(pub(crate) Rat);

impl Serialize for SerDeRat {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error> where S: serde::Serializer {
        let (numer, denom) = self.0.clone().into_signed_numerator_and_denominator();
        if denom == Natural::ONE {
            numer.to_string().serialize(serializer)
        } else {
            format!("{numer}/{denom}").serialize(serializer)
        }
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

#[derive(Clone, Debug)]
#[derive(Serialize, Deserialize)]
pub(crate) struct SerDeFold {
    #[serde(default = "crate::fold::version")]
    #[serde(with = "crate::ser_de::file_spec")]
    pub file_spec: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_creator: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_author: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_title: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_description: Option<String>,
    #[serde(default, skip_serializing_if = "IndexSet::is_empty")]
    pub file_classes: IndexSet<FileClass>,
    #[serde(default, skip_serializing_if = "TiVec::is_empty")]
    pub file_frames: TiVec<FrameIndex, Frame>,
}

#[derive(Clone, Debug)]
#[derive(Serialize, Deserialize)]
pub(crate) struct SerDeFrame {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) frame_author: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) frame_title: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) frame_description: Option<String>,
    #[serde(default, skip_serializing_if = "IndexSet::is_empty")]
    pub(crate) frame_classes: IndexSet<FrameClass>,
    #[serde(default, skip_serializing_if = "IndexSet::is_empty")]
    pub(crate) frame_attributes: IndexSet<FrameAttribute>,
    #[serde(default, skip_serializing_if = "FrameUnit::is_default")]
    pub(crate) frame_unit: FrameUnit,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(with = "crate::ser_de::basis")]
    #[serde(rename = "frame_exact:basis")]
    pub(crate) basis: Option<ArcBasis>,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "vertices_coords")]
    pub(crate) vertices_coords_f64: Option<TiVec<Vertex, Vec<f64>>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "vertices_exact:coords")]
    pub(crate) vertices_coords_exact: Option<TiVec<Vertex, Vec<Vec<SerDeRat>>>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) vertices_vertices: Option<TiVec<Vertex, Vec<Vertex>>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) vertices_edges: Option<TiVec<Vertex, Vec<Edge>>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) vertices_faces: Option<TiVec<Vertex, Vec<Option<Face>>>>,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) edges_vertices: Option<TiVec<Edge, [Vertex; 2]>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) edges_faces: Option<TiVec<Edge, Vec<Option<Face>>>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) edges_assignment: Option<TiVec<Edge, EdgeAssignment>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edges_foldAngle")]
    pub(crate) edges_fold_angle_f64: Option<TiVec<Edge, f64>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edges_exact:foldAngle")]
    pub(crate) edges_fold_angle_exact: Option<TiVec<Edge, (i32, Option<Vec<SerDeRat>>)>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edges_length")]
    pub(crate) edges_length_f64: Option<TiVec<Edge, f64>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edges_exact:length2")]
    pub(crate) edges_length2_exact: Option<TiVec<Edge, Vec<SerDeRat>>>,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) faces_vertices: Option<TiVec<Face, Vec<Vertex>>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) faces_edges: Option<TiVec<Face, Vec<Edge>>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) faces_faces: Option<TiVec<Face, Vec<Option<Face>>>>,
    
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "faceOrders")]
    pub(crate) face_orders: Option<Vec<(Face, Face, FaceOrder)>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "edgeOrders")]
    pub(crate) edge_orders: Option<Vec<(Edge, Edge, EdgeOrder)>>,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) frame_parent: Option<FrameIndex>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) frame_inherit: Option<bool>,
    #[serde(flatten)]
    pub(crate) frame_custom: Map<String, Value>,
}

impl Default for SerDeFrame {
    fn default() -> Self {
        serde_json::de::from_str("{}").unwrap()
    }
}

fn based_expr_from_de(coeffs: Vec<SerDeRat>, basis: ArcBasis) -> Result<BasedExpr, String> {
    let coeffs = coeffs.into_iter().map(|s| s.0).collect::<Vec<_>>();
    BasedExpr::try_from_coeffs_and_basis(coeffs, basis).ok_or("too many coefficients for the basis".to_owned())
}

impl TryFrom<SerDeFrame> for Frame {
    type Error = String;

    fn try_from(mut value: SerDeFrame) -> Result<Self, Self::Error> {
        if (value.vertices_coords_exact.is_some() || value.edges_fold_angle_exact.is_some() || value.edges_length2_exact.is_some())
            && value.basis.is_none()
        {
            Err("exact coordinates without a basis")?
        }
        let frame_custom = std::mem::take(&mut value.frame_custom);
        let mut vertices_custom = IndexMap::new();
        let mut edges_custom = IndexMap::new();
        let mut faces_custom = IndexMap::new();
        let mut other_custom = IndexMap::new();
        for (key, value) in frame_custom {
            if key.starts_with("vertices_") {
                let array = if let Value::Array(v) = value { v } else { Err(format!("{key} must be an array"))? };
                vertices_custom.insert(key.replacen("vertices_", "", 1), array.into());
            } else if key.starts_with("edges_") {
                let array = if let Value::Array(v) = value { v } else { Err(format!("{key} must be an array"))? };
                edges_custom.insert(key.replacen("edges_", "", 1), array.into());
            } else if key.starts_with("faces_") {
                let array = if let Value::Array(v) = value { v } else { Err(format!("{key} must be an array"))? };
                faces_custom.insert(key.replacen("faces_", "", 1), array.into());
            } else {
                other_custom.insert(key, value);
            }
        }
        let vertices_meta = value.num_vertices(&vertices_custom)
            .map_err(|e| e.into_iter().map(|e| format!("{e}\n")).collect::<String>())?;
        let edges_meta = value.num_edges(&edges_custom)
            .map_err(|e| e.into_iter().map(|e| format!("{e}\n")).collect::<String>())?;
        let faces_meta = value.num_faces(&faces_custom)
            .map_err(|e| e.into_iter().map(|e| format!("{e}\n")).collect::<String>())?;
        let num_vertices = vertices_meta.as_ref().map(|&(n, _)| n).unwrap_or(0);
        let (value, vertices_edges, edges_vertices, edges_faces, faces_half_edges) =
            match value.extract_topology(vertices_meta, edges_meta, faces_meta) {
            Ok(t) => t,
            Err(e) => Err(e.into_iter().map(|e| format!("{e}\n")).collect::<String>())?
        };

        let mut frame = Self {
            frame_author: value.frame_author,
            frame_title: value.frame_title,
            frame_description: value.frame_description,
            frame_classes: value.frame_classes,
            frame_attributes: value.frame_attributes,
            frame_unit: value.frame_unit,
            basis: value.basis.clone(),
            num_vertices,
            vertices_coords_f64: value.vertices_coords_f64
                .map(|vs| {
                    let num_rows = vs.iter().map(|v| v.len()).max().unwrap_or(0);
                    let num_cols = vs.len();
                    DMatrix::from_vec(num_rows, num_cols, vs.into_iter().flat_map(|mut v| {
                        v.resize(num_rows, 0.0);
                        v
                    }).collect::<Vec<_>>())
                }),
            vertices_coords_exact: value.vertices_coords_exact
                .map(|vs| {
                    let num_rows = vs.iter().map(|v| v.len()).max().unwrap_or(0);
                    let num_cols = vs.len();
                    let basis = value.basis.as_ref().unwrap();
                    vs.into_iter()
                        .flat_map(|mut v| {
                            v.resize(num_rows, vec![]);
                            v.into_iter()
                                .map(|c| based_expr_from_de(c, basis.clone()))
                                .collect::<Vec<_>>()
                        })
                        .collect::<Result<Vec<_>, _>>()
                        .map(|vs| DMatrix::from_vec(num_rows, num_cols, vs))
                })
                .transpose()?,
            vertices_half_edges: vertices_edges,
            //vertices_fan_boundaries: None, // will be calculated if manifold
            edges_vertices,
            edges_face_corners: edges_faces,
            edges_assignment: value.edges_assignment,
            edges_fold_angle_f64: value.edges_fold_angle_f64,
            edges_length_f64: value.edges_length_f64,
            edges_fold_angle_exact: value.edges_fold_angle_exact
                .map(|ls| ls.into_iter()
                    .map(|(floor_num_turns, tan)|
                        tan.map(|tan| based_expr_from_de(tan, value.basis.as_ref().unwrap().clone()))
                            .transpose()
                            .map(|tan| Angle::new(floor_num_turns, tan)))
                    .collect::<Result<TiVec<_, _>, _>>())
                .transpose()?,
            edges_length2_exact: value.edges_length2_exact
                .map(|ls| ls.into_iter()
                    .map(|c| based_expr_from_de(c, value.basis.as_ref().unwrap().clone()))
                    .collect::<Result<TiVec<_, _>, _>>())
                .transpose()?,
            faces_half_edges,
            face_orders: value.face_orders,
            edge_orders: value.edge_orders,
            frame_parent: value.frame_parent,
            frame_inherit: value.frame_inherit,
            vertices_custom,
            edges_custom,
            faces_custom,
            other_custom,
        };

        frame = frame.check()
            .map_err(|e| e.into_iter().map(|e| format!("{e}\n")).collect::<String>())?;
        
        Ok(frame)
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
        let (value,
            vertices_vertices, vertices_edges, vertices_faces,
               edges_vertices,                    edges_faces,
               faces_vertices,    faces_edges,    faces_faces) = match value.extract_topology()
        {
            Ok(t) => t,
            Err(e) => todo!("{}", e),
        };

        Self {
            frame_author: value.frame_author,
            frame_title: value.frame_title,
            frame_description: value.frame_description,
            frame_classes: value.frame_classes,
            frame_attributes: value.frame_attributes,
            frame_unit: value.frame_unit,
            basis: value.basis,
            vertices_coords_f64: value.vertices_coords_f64
                .map(|vs| vs.column_iter()
                    .map(|v| v.into_iter().copied().collect::<Vec<_>>())
                    .collect::<TiVec<_, _>>()),
            vertices_coords_exact: value.vertices_coords_exact
                .map(|mut vs| vs.column_iter_mut()
                    .map(|v| v.into_iter().map(std::mem::take).map(based_expr_to_ser).collect::<Vec<_>>())
                    .collect::<TiVec<_, _>>()),
            vertices_vertices,
            vertices_edges,
            vertices_faces,
            edges_vertices,
            edges_faces,
            edges_assignment: value.edges_assignment,
            edges_fold_angle_f64: value.edges_fold_angle_f64,
            edges_length_f64: value.edges_length_f64,
            edges_fold_angle_exact: value.edges_fold_angle_exact
                .map(|ls| ls.into_iter()
                    .map(|angle| angle.into_turn_value_tan())
                    .map(|(turn_value, tan)|
                        (turn_value, tan.map(based_expr_to_ser)))
                    .collect::<TiVec<_, _>>()),
            edges_length2_exact: value.edges_length2_exact
                .map(|ls| ls.into_iter().map(based_expr_to_ser).collect::<TiVec<_, _>>()),
            faces_vertices,
            faces_edges,
            faces_faces,
            face_orders: value.face_orders,
            edge_orders: value.edge_orders,
            frame_parent: value.frame_parent,
            frame_inherit: value.frame_inherit,
            frame_custom: value.vertices_custom.into_iter().map(|(k, v)| (format!("vertices_{k}"), Value::Array(v.into())))
                .chain(value.edges_custom.into_iter().map(|(k, v)| (format!("edges_{k}"), Value::Array(v.into()))))
                .chain(value.faces_custom.into_iter().map(|(k, v)| (format!("faces_{k}"), Value::Array(v.into()))))
                .chain(value.other_custom)
                .collect::<Map<_, _>>()
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
    
    pub fn to_json_value(&self) -> Result<Value, serde_json::error::Error> {
        // Done like this to avoid unnecessary clones
        // (instead of converting to SerDeFold and serializing that)
        let spec = self.file_spec.parse::<f64>()
            .map_err(|e| serde_json::error::Error::custom(format!("unfortunately spec must be convertible to a number: {e:?}")))?;
        let mut ser_de_fold = vec![
                                                      Some(("file_spec"       .to_owned(), serde_json::to_value(&spec)?)),
            if let Some(c) = &self.file_creator     { Some(("file_creator"    .to_owned(), serde_json::to_value(c)?)) } else { None },
            if let Some(c) = &self.file_author      { Some(("file_author"     .to_owned(), serde_json::to_value(c)?)) } else { None },
            if let Some(c) = &self.file_title       { Some(("file_title"      .to_owned(), serde_json::to_value(c)?)) } else { None },
            if let Some(c) = &self.file_description { Some(("file_description".to_owned(), serde_json::to_value(c)?)) } else { None },
            if self.file_classes.len() > 0          { Some(("file_classes"    .to_owned(), serde_json::to_value(&self.file_classes)?)) } else { None },
            if self.file_frames .len() > 1          { Some(("file_frames"     .to_owned(), serde_json::to_value(&self.file_frames[FrameIndex(1)..])?)) } else { None },
        ].into_iter().flatten().collect::<Vec<_>>();
        let key_frame = serde_json::to_value(&self.file_frames.first().unwrap())?;
        let key_frame = if let Value::Object(k) = key_frame { k } else { unreachable!() };
        let file_custom = serde_json::to_value(&self.file_custom)?;
        let file_custom = if let Value::Object(k) = file_custom { k } else { unreachable!() };
        ser_de_fold.extend(key_frame);
        ser_de_fold.extend(file_custom);
        let ser_de_fold = ser_de_fold.into_iter().collect::<Map<_, _>>();
        Ok(Value::Object(ser_de_fold))
    }
    
    pub fn to_json_string(&self) -> Result<String, serde_json::error::Error> {
        pretty::to_string_fold(&self.to_json_value()?)
    }

    /// #[serde(flatten)] is glitched for structs, so we have to work around it
    pub fn from_str(s: &str) -> Result<Self, serde_json::error::Error> {
        // with a double deserialize!
        let mut fold = serde_json::de::from_str::<SerDeFold>(s)?;
        let key_frame = serde_json::de::from_str::<Frame>(s)?;
        fold.file_frames.insert(FrameIndex(0), key_frame);
        let mut fold = Fold {
            file_spec: fold.file_spec,
            file_creator: fold.file_creator,
            file_author: fold.file_author,
            file_title: fold.file_title,
            file_description: fold.file_description,
            file_classes: fold.file_classes,
            file_frames: fold.file_frames,
            file_custom: IndexMap::new()
        };

        // Move custom fields appropriately
        let mut file_custom = IndexMap::new();
        let mut frame_custom = IndexMap::new();
        for keys in std::mem::take(&mut fold.key_frame_mut().other_custom) {
            if keys.0.starts_with("frame_") {
                frame_custom.insert(keys.0, keys.1);
            } else if keys.0.contains(":") {
                file_custom.insert(keys.0, keys.1);
            }
        }
        
        fold.file_custom = file_custom;
        fold.key_frame_mut().other_custom = frame_custom;
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
    use indexmap::{indexset, IndexMap, IndexSet};
    use exact_number::{based_expr, basis::{ArcBasis, Basis}, sqrt_expr, Angle};
    use nalgebra::DMatrix;
    use serde_json::{json, Value};
    use typed_index_collections::ti_vec;

    use crate::{fold::{EdgeAssignment, FileClass, Fold, Frame, FrameAttribute, FrameClass, FrameIndex, FrameUnit}};
    use crate::fold::{Vertex as V, Edge as E, HalfEdge as H, Face as F, FaceCorner as C};

    #[test]
    fn test_load_partial_metadata() {
        let fold = Fold::from_file("test/x.fold").unwrap();
        assert_eq!(fold.file_spec, "1.1");
        assert_eq!(fold.file_creator, Some("oriedita".to_owned()));
        assert_eq!(fold.file_author, None);
        assert_eq!(fold.file_title, None);
        assert_eq!(fold.file_description, None);
        assert_eq!(fold.file_classes, IndexSet::new());
        assert_eq!(fold.file_custom, vec![("oriedita:version".to_owned(), Value::String("1.1.3".to_owned()))].into_iter().collect::<IndexMap<_, _>>());
    }

    #[test]
    fn test_load_full_metadata() {
        let fold = Fold::from_file("test/x-all-file-metadata.fold").unwrap();
        assert_eq!(fold.file_spec, "1.1");
        assert_eq!(fold.file_creator, Some("oriedita".to_owned()));
        assert_eq!(fold.file_author, Some("hayastl".to_owned()));
        assert_eq!(fold.file_title, Some("x".to_owned()));
        assert_eq!(fold.file_description, Some("a simple cross pattern".to_owned()));
        assert_eq!(fold.file_classes, indexset![FileClass::SingleModel]);
        assert_eq!(fold.file_custom, vec![("oriedita:version".to_owned(), Value::String("1.1.3".to_owned()))].into_iter().collect::<IndexMap<_, _>>());
    }

    #[test]
    fn test_load_exact_coords() {
        let fold = Fold::from_file("test/x-exact-coords.fold").unwrap();
        let basis = Basis::new_arc(vec![sqrt_expr!(1)]);
        let coords = DMatrix::from_vec(2, 5, vec![
            based_expr!(-200), based_expr!(-200),
            based_expr!(-200), based_expr!(200),
            based_expr!(200), based_expr!(-200),
            based_expr!(200), based_expr!(200),
            based_expr!(0), based_expr!(0),
        ]);
        assert_eq!(fold.key_frame().basis, Some(basis));
        assert_eq!(fold.key_frame().vertices_coords_exact, Some(coords));
    }

    #[test]
    fn test_load_exact_coords_bird_base() {
        let fold = Fold::from_file("test/x-exact-coords-bird-base.fold").unwrap();
        let basis = Basis::new_arc(vec![sqrt_expr!(1), sqrt_expr!(2)]);
        let coords = DMatrix::from_vec(2, 9, vec![
            based_expr!(0 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2),
            based_expr!(1 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2),
            based_expr!(1 + 0 sqrt 2), based_expr!(1 + 0 sqrt 2),
            based_expr!(0 + 0 sqrt 2), based_expr!(1 + 0 sqrt 2),
            based_expr!(1/2 + 0 sqrt 2), based_expr!(1/2 + 0 sqrt 2),
            based_expr!(1/2 + 0 sqrt 2), based_expr!(-1/2 + 1/2 sqrt 2),
            based_expr!(3/2 - 1/2 sqrt 2), based_expr!(1/2 + 0 sqrt 2),
            based_expr!(1/2 + 0 sqrt 2), based_expr!(3/2 - 1/2 sqrt 2),
            based_expr!(-1/2 + 1/2 sqrt 2), based_expr!(1/2 + 0 sqrt 2),
        ]);
        assert_eq!(fold.key_frame().basis, Some(basis));
        assert_eq!(fold.key_frame().vertices_coords_exact, Some(coords));
    }

    #[test]
    fn test_load_all_fields() {
        let fold = Fold::from_file("test/x-all-file-data.fold").unwrap();
        
        assert_eq!(fold.file_spec, "1.1");
        assert_eq!(fold.file_creator, Some("oriedita".to_owned()));
        assert_eq!(fold.file_author, Some("hayastl".to_owned()));
        assert_eq!(fold.file_title, Some("x".to_owned()));
        assert_eq!(fold.file_description, Some("a simple cross pattern".to_owned()));
        assert_eq!(fold.file_classes, indexset![FileClass::SingleModel]);
        assert_eq!(fold.file_frames[FrameIndex(1)].frame_title, Some("Nothing".to_owned()));
        assert_eq!(fold.file_frames[FrameIndex(1)].frame_parent, Some(FrameIndex(0)));
        assert_eq!(fold.file_custom, vec![("oriedita:version".to_owned(), Value::String("1.1.3".to_owned()))].into_iter().collect::<IndexMap<_, _>>());
        assert_eq!(fold.key_frame().frame_author, Some("hayastl".to_owned()));
        assert_eq!(fold.key_frame().frame_title, Some("The Key Frame".to_owned()));
        assert_eq!(fold.key_frame().frame_description, Some("the crease pattern".to_owned()));
        assert_eq!(fold.key_frame().frame_classes, indexset![FrameClass::CreasePattern]);
        assert_eq!(fold.key_frame().frame_attributes, indexset![FrameAttribute::_2D, FrameAttribute::Manifold, FrameAttribute::Orientable,
            FrameAttribute::NonSelfTouching, FrameAttribute::NonSelfIntersecting,
            FrameAttribute::NoCuts, FrameAttribute::NoJoins, FrameAttribute::ConvexFaces]);
        assert_eq!(fold.key_frame().frame_unit, FrameUnit::Millimeter);
        assert_eq!(fold.key_frame().basis, Some(Basis::new_arc(vec![sqrt_expr!(1)])));
        assert_eq!(fold.key_frame().vertices_coords_f64, Some(DMatrix::from_vec(2, 5, vec![
            -200.0, -200.0,
            -200.0, 200.0,
            200.0, 200.0,
            200.0, -200.0,
            0.0, 0.0,
        ])));
        assert_eq!(fold.key_frame().vertices_coords_exact, Some(DMatrix::from_vec(2, 5, vec![
            based_expr!(-200), based_expr!(-200),
            based_expr!(-200), based_expr!(200),
            based_expr!(200), based_expr!(200),
            based_expr!(200), based_expr!(-200),
            based_expr!(0), based_expr!(0),
        ])));
        assert_eq!(fold.key_frame().vertices_half_edges, Some(ti_vec![
            vec![H(0), H(8), H(7)],
            vec![H(2), H(10), H(1)],
            vec![H(4), H(12), H(3)],
            vec![H(6), H(14), H(5)],
            vec![H(9), H(11), H(13), H(15)],
        ]));
        assert_eq!(fold.key_frame().edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(4)],
            [V(1), V(4)],
            [V(2), V(4)],
            [V(3), V(4)],
        ]));
        assert_eq!(fold.key_frame().edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(1), 0)], vec![]],
            [vec![C(F(2), 0)], vec![]],
            [vec![C(F(3), 0)], vec![]],
            [vec![C(F(3), 1)], vec![C(F(0), 2)]],
            [vec![C(F(0), 1)], vec![C(F(1), 2)]],
            [vec![C(F(1), 1)], vec![C(F(2), 2)]],
            [vec![C(F(2), 1)], vec![C(F(3), 2)]],
        ]));
        assert_eq!(fold.key_frame().edges_assignment, Some(ti_vec![
            EdgeAssignment::Boundary,
            EdgeAssignment::Boundary,
            EdgeAssignment::Boundary,
            EdgeAssignment::Boundary,
            EdgeAssignment::Mountain,
            EdgeAssignment::Mountain,
            EdgeAssignment::Mountain,
            EdgeAssignment::Valley,
        ]));
        assert_eq!(fold.key_frame().edges_fold_angle_f64, Some(ti_vec![
            0.0,
            0.0,
            0.0,
            0.0,
            -180.0,
            -180.0,
            -180.0,
            180.0,
        ]));
        assert_eq!(fold.key_frame().edges_fold_angle_exact, Some(ti_vec![
            Angle::new( 0, Some(based_expr!(0))),
            Angle::new( 0, Some(based_expr!(0))),
            Angle::new( 0, Some(based_expr!(0))),
            Angle::new( 0, Some(based_expr!(0))),
            Angle::new(-1, Some(based_expr!(0))),
            Angle::new(-1, Some(based_expr!(0))),
            Angle::new(-1, Some(based_expr!(0))),
            Angle::new( 1, Some(based_expr!(0))),
        ]));
        assert_eq!(fold.key_frame().edges_length_f64, Some(ti_vec![
            200.0,
            200.0,
            200.0,
            200.0,
            141.4213562373095,
            141.4213562373095,
            141.4213562373095,
            141.4213562373095,
        ]));
        assert_eq!(fold.key_frame().edges_length2_exact, Some(ti_vec![
            based_expr!(40000),
            based_expr!(40000),
            based_expr!(40000),
            based_expr!(40000),
            based_expr!(20000),
            based_expr!(20000),
            based_expr!(20000),
            based_expr!(20000),
        ]));
        assert_eq!(fold.key_frame().faces_half_edges, Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ]));
        assert_eq!(fold.key_frame().face_orders, Some(vec![]));
        assert_eq!(fold.key_frame().edge_orders, Some(vec![]));
        assert_eq!(fold.key_frame().frame_parent, None);
        assert_eq!(fold.key_frame().frame_inherit, Some(true));
        assert_eq!(fold.key_frame().vertices_custom, vec![
            ("test:degree".to_owned(), ti_vec![
                json!(3),
                json!(3),
                json!(3),
                json!(3),
                json!(4),
            ])
        ].into_iter().collect::<IndexMap<_, _>>());
        assert_eq!(fold.key_frame().edges_custom, vec![
            ("test:color".to_owned(), ti_vec![
                json!("red"),
                json!("green"),
                json!("red"),
                json!("green"),
                json!("blue"),
                json!("yellow"),
                json!("blue"),
                json!("yellow"),
            ])
        ].into_iter().collect::<IndexMap<_, _>>());
        assert_eq!(fold.key_frame().faces_custom, vec![
            ("test:circumcenter".to_owned(), ti_vec![
                json!([100, 0]),
                json!([200, 100]),
                json!([100, 200]),
                json!([0, 100]),
            ])
        ].into_iter().collect::<IndexMap<_, _>>());
        assert_eq!(fold.key_frame().other_custom, vec![("frame_test:foo".to_owned(), Value::String("bar".to_owned()))].into_iter().collect::<IndexMap<_, _>>());
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

    #[test]
    fn test_serialize_simple() {
        let fold = Fold {
            file_spec: "1.2".to_owned(),
            file_creator: Some("rust".to_owned()),
            file_frames: ti_vec![
                    Frame {
                    vertices_coords_f64: Some(DMatrix::from_vec(2, 4, vec![
                        0.0, 0.0,
                        1.0, 0.0,
                        1.0, 1.0,
                        0.0, 1.0,
                    ])),
                    ..Default::default()
                },
            ],
            ..Default::default()
        };
        let expected = json!({
            "file_spec": 1.2,
            "file_creator": "rust",
            "vertices_coords": [
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [0.0, 1.0]
            ]
        });
        let result = fold.to_json_value().unwrap();
        assert_eq!(result, expected);
        println!("test_serialize_simple output: {}", fold.to_json_string().unwrap());
    }

    #[test]
    fn test_serialize_exact_coords() {
        let fold = Fold {
            file_spec: "1.2".to_owned(),
            file_creator: Some("rust".to_owned()),
            file_frames: ti_vec![
                Frame {
                    basis: Some(Basis::new_arc(vec![sqrt_expr!(1), sqrt_expr!(2)])),
                    vertices_coords_exact: Some(DMatrix::from_vec(2, 9, vec![
                        based_expr!(0 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2),
                        based_expr!(1 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2),
                        based_expr!(1 + 0 sqrt 2), based_expr!(1 + 0 sqrt 2),
                        based_expr!(0 + 0 sqrt 2), based_expr!(1 + 0 sqrt 2),
                        based_expr!(1/2 + 0 sqrt 2), based_expr!(1/2 + 0 sqrt 2),
                        based_expr!(1/2 + 0 sqrt 2), based_expr!(-1/2 + 1/2 sqrt 2),
                        based_expr!(3/2 - 1/2 sqrt 2), based_expr!(1/2 + 0 sqrt 2),
                        based_expr!(1/2 + 0 sqrt 2), based_expr!(3/2 - 1/2 sqrt 2),
                        based_expr!(-1/2 + 1/2 sqrt 2), based_expr!(1/2 + 0 sqrt 2),
                    ])),
                    ..Default::default()
                },
            ],
            ..Default::default()
        };
        let expected = json!({
            "file_spec": 1.2,
            "file_creator": "rust",
            "frame_exact:basis": ["√1", "√2"],
            "vertices_exact:coords": [
                [["0", "0"], ["0", "0"]],
                [["1", "0"], ["0", "0"]],
                [["1", "0"], ["1", "0"]],
                [["0", "0"], ["1", "0"]],
                [["1/2", "0"], ["1/2", "0"]],
                [["1/2", "0"], ["-1/2", "1/2"]],
                [["3/2", "-1/2"], ["1/2", "0"]],
                [["1/2", "0"], ["3/2", "-1/2"]],
                [["-1/2", "1/2"], ["1/2", "0"]],
            ]
        });
        let result = fold.to_json_value().unwrap();
        assert_eq!(result, expected);
        println!("test_serialize_exact_coords output: {}", fold.to_json_string().unwrap());
    }

    #[test]
    fn test_serialize_multiple_frames() {
        let fold = Fold {
            file_spec: "1.2".to_owned(),
            file_creator: Some("rust".to_owned()),
            file_frames: ti_vec![
                Frame {
                    vertices_coords_f64: Some(DMatrix::from_vec(2, 4, vec![
                        0.0, 0.0,
                        1.0, 0.0,
                        1.0, 1.0,
                        0.0, 1.0,
                    ])),
                    ..Default::default()
                },
                Frame {
                    frame_title: Some("Scaled by 2".to_owned()),
                    vertices_coords_f64: Some(DMatrix::from_vec(2, 4, vec![
                        0.0, 0.0,
                        2.0, 0.0,
                        2.0, 2.0,
                        0.0, 2.0,
                    ])),
                    ..Default::default()
                },
                Frame {
                    frame_title: Some("Scaled by 3".to_owned()),
                    vertices_coords_f64: Some(DMatrix::from_vec(2, 4, vec![
                        0.0, 0.0,
                        3.0, 0.0,
                        3.0, 3.0,
                        0.0, 3.0,
                    ])),
                    ..Default::default()
                }
            ],
            ..Default::default()
        };
        let string = fold.to_json_string().unwrap();
        let num_lines = string.lines().count();
        assert_eq!(num_lines, 30); // make sure only the appropriate arrays are inlined
        println!("test_serialize_multiple_frames output: {}", string);
    }

    #[test]
    fn test_serialize_all_fields() {
        // Let's just make sure a round trip preserves equality.
        let fold = Fold {
            file_spec: "1.1".to_owned(),
            file_creator: Some("oriedita".to_owned()),
            file_author: Some("hayastl".to_owned()),
            file_title: Some("x".to_owned()),
            file_description: Some("a simple cross pattern".to_owned()),
            file_classes: indexset![FileClass::SingleModel],
            file_frames: ti_vec![
                Frame {
                    frame_author: Some("hayastl".to_owned()),
                    frame_title: Some("The Key Frame".to_owned()),
                    frame_description: Some("the crease pattern".to_owned()),
                    frame_classes: indexset![FrameClass::CreasePattern],
                    frame_attributes: indexset![FrameAttribute::_2D, FrameAttribute::Manifold, FrameAttribute::Orientable,
                        FrameAttribute::NonSelfTouching, FrameAttribute::NonSelfIntersecting,
                        FrameAttribute::NoCuts, FrameAttribute::NoJoins, FrameAttribute::ConvexFaces],
                    frame_unit: FrameUnit::Millimeter,
                    basis: Some(Basis::new_arc(vec![sqrt_expr!(1)])),
                    num_vertices: 5,
                    vertices_coords_f64: Some(DMatrix::from_vec(2, 5, vec![
                        -200.0, -200.0,
                        -200.0, 200.0,
                        200.0, 200.0,
                        200.0, -200.0,
                        0.0, 0.0,
                    ])),
                    vertices_coords_exact: Some(DMatrix::from_vec(2, 5, vec![
                        based_expr!(-200), based_expr!(-200),
                        based_expr!(-200), based_expr!(200),
                        based_expr!(200), based_expr!(200),
                        based_expr!(200), based_expr!(-200),
                        based_expr!(0), based_expr!(0),
                    ])),
                    vertices_half_edges: Some(ti_vec![
                        vec![H(0), H(8), H(7)],
                        vec![H(2), H(10), H(1)],
                        vec![H(4), H(12), H(3)],
                        vec![H(6), H(14), H(5)],
                        vec![H(9), H(11), H(13), H(15)],
                    ]),
                    //vertices_fan_boundaries: Some(ti_vec![
                    //    vec![0, 3],
                    //    vec![0, 3],
                    //    vec![0, 3],
                    //    vec![0, 3],
                    //    vec![]
                    //]),
                    edges_vertices: Some(ti_vec![
                        [V(0), V(1)],
                        [V(1), V(2)],
                        [V(2), V(3)],
                        [V(3), V(0)],
                        [V(0), V(4)],
                        [V(1), V(4)],
                        [V(2), V(4)],
                        [V(3), V(4)],
                    ]),
                    edges_face_corners: Some(ti_vec![
                        [vec![C(F(0), 0)], vec![]],
                        [vec![C(F(1), 0)], vec![]],
                        [vec![C(F(2), 0)], vec![]],
                        [vec![C(F(3), 0)], vec![]],
                        [vec![C(F(3), 1)], vec![C(F(0), 2)]],
                        [vec![C(F(0), 1)], vec![C(F(1), 2)]],
                        [vec![C(F(1), 1)], vec![C(F(2), 2)]],
                        [vec![C(F(2), 1)], vec![C(F(3), 2)]],
                    ]),
                    edges_assignment: Some(ti_vec![
                        EdgeAssignment::Boundary,
                        EdgeAssignment::Boundary,
                        EdgeAssignment::Boundary,
                        EdgeAssignment::Boundary,
                        EdgeAssignment::Mountain,
                        EdgeAssignment::Mountain,
                        EdgeAssignment::Mountain,
                        EdgeAssignment::Valley,
                    ]),
                    edges_fold_angle_f64: Some(ti_vec![
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        -180.0,
                        -180.0,
                        -180.0,
                        180.0,
                    ]),
                    edges_fold_angle_exact: Some(ti_vec![
                        Angle::new( 0, Some(based_expr!(0))),
                        Angle::new( 0, Some(based_expr!(0))),
                        Angle::new( 0, Some(based_expr!(0))),
                        Angle::new( 0, Some(based_expr!(0))),
                        Angle::new(-1, Some(based_expr!(0))),
                        Angle::new(-1, Some(based_expr!(0))),
                        Angle::new(-1, Some(based_expr!(0))),
                        Angle::new( 1, Some(based_expr!(0))),
                    ]),
                    edges_length_f64: Some(ti_vec![
                        200.0,
                        200.0,
                        200.0,
                        200.0,
                        141.4213562373095,
                        141.4213562373095,
                        141.4213562373095,
                        141.4213562373095,
                    ]),
                    edges_length2_exact: Some(ti_vec![
                        based_expr!(40000),
                        based_expr!(40000),
                        based_expr!(40000),
                        based_expr!(40000),
                        based_expr!(20000),
                        based_expr!(20000),
                        based_expr!(20000),
                        based_expr!(20000),
                    ]),
                    faces_half_edges: Some(ti_vec![
                        vec![H(0), H(10), H(9)],
                        vec![H(2), H(12), H(11)],
                        vec![H(4), H(14), H(13)],
                        vec![H(6), H(8), H(15)],
                    ]),
                    face_orders: Some(vec![]),
                    edge_orders: Some(vec![]),
                    frame_parent: None,
                    frame_inherit: Some(true),
                    vertices_custom: vec![
                        ("test:degree".to_owned(), ti_vec![
                            json!(3),
                            json!(3),
                            json!(3),
                            json!(3),
                            json!(4),
                        ])
                    ].into_iter().collect::<IndexMap<_, _>>(),
                    edges_custom: vec![
                        ("test:color".to_owned(), ti_vec![
                            json!("red"),
                            json!("green"),
                            json!("red"),
                            json!("green"),
                            json!("blue"),
                            json!("yellow"),
                            json!("blue"),
                            json!("yellow"),
                        ])
                    ].into_iter().collect::<IndexMap<_, _>>(),
                    faces_custom: vec![
                        ("test:circumcenter".to_owned(), ti_vec![
                            json!([100, 0]),
                            json!([200, 100]),
                            json!([100, 200]),
                            json!([0, 100]),
                        ])
                    ].into_iter().collect::<IndexMap<_, _>>(),
                    other_custom: vec![("frame_test:foo".to_owned(), Value::String("bar".to_owned()))].into_iter().collect::<IndexMap<_, _>>(),
                },
                Frame {
                    frame_title: Some("Nothing".to_owned()),
                    frame_parent: Some(FrameIndex(0)),
                    ..Default::default()
                }
            ],
            file_custom: vec![("oriedita:version".to_owned(), Value::String("1.1.3".to_owned()))].into_iter().collect::<IndexMap<_, _>>(),
        };
        let string = fold.to_json_string().unwrap();
        println!("test_serialize_all_fields output: {}", string);
        //let round_trip_fold = Fold::from_str(&string).unwrap();
        //assert_eq!(fold, round_trip_fold);
    }

}