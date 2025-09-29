//use std::{collections::HashMap, fs::File, path::Path};
//
//use json::{object::Object, Array, JsonValue};
//
//use crate::fold::{FileClass, Fold, Frame, FrameAttribute, FrameClass, FrameUnit};
//
//#[derive(Debug)]
//pub enum FoldFileError {
//    IoError(std::io::Error),
//    JsonError(json::Error),
//    InvalidString(String, JsonValue),
//    InvalidFileClass(String),
//    InvalidFrameClass(String),
//    InvalidFrameAttribute(String),
//    InvalidFrameUnit(String),
//}
//
//impl FileClass {
//    fn from_str(str: String) -> Result<Self, FoldFileError> {
//        if &str == "singleModel" { Ok(Self::SingleModel) } else
//        if &str == "multiModel" { Ok(Self::MultiModel) } else
//        if &str == "animation" { Ok(Self::Animation) } else
//        if &str == "diagrams" { Ok(Self::Diagrams) } else
//        if str.contains(":") { Ok(Self::Custom(str)) } else
//        { Err(FoldFileError::InvalidFileClass(str)) }
//    }
//
//    fn from_json<'a>(val: impl Iterator<Item = &'a JsonValue>) -> Result<Vec<Self>, FoldFileError> {
//        let classes = val
//            .map(|v| Self::from_str(v.to_string()))
//            .collect::<Result<Vec<_>, _>>();
//        classes
//    }
//}
//
//impl FrameClass {
//    fn from_str(str: String) -> Result<Self, FoldFileError> {
//        if &str == "creasePattern" { Ok(Self::CreasePattern) } else
//        if &str == "foldedForm" { Ok(Self::FoldedForm) } else
//        if &str == "graph" { Ok(Self::Graph) } else
//        if &str == "linkage" { Ok(Self::Linkage) } else
//        if str.contains(":") { Ok(Self::Custom(str)) } else
//        { Err(FoldFileError::InvalidFrameClass(str)) }
//    }
//
//    fn from_json<'a>(val: impl Iterator<Item = &'a JsonValue>) -> Result<Vec<Self>, FoldFileError> {
//        let classes = val
//            .map(|v| Self::from_str(v.to_string()))
//            .collect::<Result<Vec<_>, _>>();
//        classes
//    }
//}
//
//impl FrameAttribute {
//    fn from_str(str: String) -> Result<Self, FoldFileError> {
//        if &str == "2D" { Ok(Self::_2D) } else
//        if &str == "3D" { Ok(Self::_3D) } else
//        if &str == "abstract" { Ok(Self::Abstract) } else
//        if &str == "manifold" { Ok(Self::Manifold(true)) } else
//        if &str == "nonManifold" { Ok(Self::Manifold(false)) } else
//        if &str == "orientable" { Ok(Self::Orientable(true)) } else
//        if &str == "nonOrientable" { Ok(Self::Orientable(false)) } else
//        if &str == "selfTouching" { Ok(Self::SelfTouching(true)) } else
//        if &str == "nonSelfTouching" { Ok(Self::SelfTouching(false)) } else
//        if &str == "selfIntersecting" { Ok(Self::SelfIntersecting(true)) } else
//        if &str == "nonSelfIntersecting" { Ok(Self::SelfIntersecting(false)) } else
//        if &str == "cuts" { Ok(Self::Cuts(true)) } else
//        if &str == "nonCuts" { Ok(Self::Cuts(false)) } else
//        if &str == "joins" { Ok(Self::Joins(true)) } else
//        if &str == "nonJoins" { Ok(Self::Joins(false)) } else
//        if &str == "convexFaces" { Ok(Self::ConvexFaces(true)) } else
//        if &str == "nonConvexFaces" { Ok(Self::ConvexFaces(false)) } else
//        if str.contains(":") { Ok(Self::Custom(str)) } else
//        { Err(FoldFileError::InvalidFrameAttribute(str)) }
//    }
//
//    fn from_json<'a>(val: impl Iterator<Item = &'a JsonValue>) -> Result<Vec<Self>, FoldFileError> {
//        let attributes = val
//            .map(|v| Self::from_str(v.to_string()))
//            .collect::<Result<Vec<_>, _>>();
//        attributes
//    }
//}
//
//impl FrameUnit {
//    fn from_str(str: String) -> Result<Self, FoldFileError> {
//        if &str == "unit" { Ok(Self::Unit) } else
//        if &str == "in" { Ok(Self::Inch) } else
//        if &str == "pt" { Ok(Self::Point) } else
//        if &str == "m" { Ok(Self::Meter) } else
//        if &str == "cm" { Ok(Self::Centimeter) } else
//        if &str == "mm" { Ok(Self::Millimeter) } else
//        if &str == "um" { Ok(Self::Micrometer) } else
//        if &str == "nm" { Ok(Self::Nanometer) } else
//        { Err(FoldFileError::InvalidFrameUnit(str)) }
//    }
//}
//
//impl<N> Frame<N> {
//    fn from_obj(obj: Object) -> Result<Self, FoldFileError> {
//        let custom = obj.iter()
//            .filter(|(k, _)| !k.starts_with("file_") && k.contains(":"))
//            .map(|(k, v)| (k.to_owned(), v.to_string()))
//            .collect::<HashMap<_, _>>();
//
//        let frame = Frame {
//            frame_author: obj.get("frame_author").map(|v| v.to_string()),
//            frame_title: obj.get("frame_title").map(|v| v.to_string()),
//            frame_description: obj.get("frame_description").map(|v| v.to_string()),
//            frame_classes: obj.get("frame_classes").map_or(Ok(vec![]), |v| FrameClass::from_json(v.members()))?,
//            frame_attributes: obj.get("frame_attributes").map_or(Ok(vec![]), |v| FrameAttribute::from_json(v.members()))?,
//            frame_unit: obj.get("frame_unit").map_or(Ok(FrameUnit::Unit), |v| FrameUnit::from_str(v.to_string()))?,
//            vertices_coords: None,
//            vertices_vertices: None,
//            vertices_edges: None,
//            vertices_faces: None,
//            edges_vertices: None,
//            edges_faces: None,
//            edges_assignment: None,
//            edges_fold_angle: None,
//            edges_exact_fold_angle: None,
//            edges_length: None,
//            faces_vertices: None,
//            faces_edges: None,
//            faces_faces: None,
//            face_orders: None,
//            edge_orders: None,
//            frame_parent: obj.get("frame_parent").and_then(|v| v.as_usize()),
//            frame_inherit: obj.get("frame_inherit").and_then(|v| v.as_bool()),
//            frame_custom: custom,
//        };
//        Ok(frame)
//    }
//}
//
//impl<N> Fold<N> {
//    pub fn from_json(json: JsonValue) -> Result<Self, FoldFileError> {
//        let obj = Object::from_iter(json.entries().map(|(k, v)| (k, v.to_owned())));
//        let custom = obj.iter()
//            .filter(|(k, _)| k.starts_with("file_") && k.contains(":"))
//            .map(|(k, v)| (k.to_owned(), v.to_string()))
//            .collect::<HashMap<_, _>>();
//
//        let mut fold = Fold {
//            file_spec: obj.get("file_spec").map_or("1.2".to_owned(), |v| v.to_string()),
//            file_creator: obj.get("file_creator").map(|v| v.to_string()),
//            file_author: obj.get("file_author").map(|v| v.to_string()),
//            file_title: obj.get("file_title").map(|v| v.to_string()),
//            file_description: obj.get("file_description").map(|v| v.to_string()),
//            file_classes: obj.get("file_classes").map_or(Ok(vec![]), |v| FileClass::from_json(v.members()))?,
//            file_frames: vec![],
//            file_custom: custom,
//        };
//        fold.file_frames.push(Frame::from_obj(obj)?);
//        Ok(fold)
//    }
//
//    pub fn from_str(s: &str) -> Result<Self, FoldFileError> {
//        let json = json::parse(s).map_err(|e| FoldFileError::JsonError(e))?;
//        Self::from_json(json)
//    }
//
//    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, FoldFileError> {
//        let text = std::fs::read_to_string(path).map_err(|e| FoldFileError::IoError(e))?;
//        Self::from_str(&text)
//    }
//}
//
//#[cfg(test)]
//mod test {
//    use std::collections::HashMap;
//
//    use crate::fold::{FileClass, FoldF64};
//
//    #[test]
//    fn test_load_partial_metadata() {
//        let fold = FoldF64::from_file("test/x.fold").unwrap();
//        assert_eq!(fold.file_spec, "1.1");
//        assert_eq!(fold.file_creator, Some("oriedita".to_owned()));
//        assert_eq!(fold.file_author, None);
//        assert_eq!(fold.file_title, None);
//        assert_eq!(fold.file_description, None);
//        assert_eq!(fold.file_classes, vec![]);
//        assert_eq!(fold.file_custom, vec![("oriedita:version".to_owned(), "1.1.3".to_owned())].into_iter().collect::<HashMap<_, _>>());
//    }
//
//    #[test]
//    fn test_load_full_metadata() {
//        let fold = FoldF64::from_file("test/x-all-metadata.fold").unwrap();
//        assert_eq!(fold.file_spec, "1.1");
//        assert_eq!(fold.file_creator, Some("oriedita".to_owned()));
//        assert_eq!(fold.file_author, Some("hayastl".to_owned()));
//        assert_eq!(fold.file_title, Some("x".to_owned()));
//        assert_eq!(fold.file_description, Some("a simple cross pattern".to_owned()));
//        assert_eq!(fold.file_classes, vec![FileClass::SingleModel]);
//        assert_eq!(fold.file_custom, vec![("oriedita:version".to_owned(), "1.1.3".to_owned())].into_iter().collect::<HashMap<_, _>>());
//    }
//
//    #[test]
//    fn test_file_classes() {
//        let result = FileClass::from_str("str")
//    }
//}