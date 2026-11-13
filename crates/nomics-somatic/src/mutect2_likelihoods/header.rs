//! htsjdk 4.2.0 SamFileHeaderMerger ID translations, including PG ancestry.
//! LinkedHashMap iteration is represented by explicit insertion-ordered vectors.
use nomics_io::bam::Header;
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    io,
};
type Id = Vec<u8>;
type RecordGroups = Vec<(Id, Vec<(BTreeMap<Id, Id>, Vec<usize>)>)>;
type Translations = Vec<HashMap<Id, Id>>;
#[derive(Clone)]
struct Entry {
    file: usize,
    id: Id,
    attributes: BTreeMap<Id, Id>,
}
fn base36(mut n: usize) -> String {
    let mut bytes = Vec::new();
    loop {
        bytes.push(b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[n % 36]);
        n /= 36;
        if n == 0 {
            break;
        }
    }
    bytes.reverse();
    String::from_utf8(bytes).unwrap()
}
pub(super) fn translations(headers: &[Header], kind: &[u8]) -> io::Result<Translations> {
    let mut pending = Vec::new();
    for (file, h) in headers.iter().enumerate() {
        let mut seen = HashSet::new();
        for line in h
            .text
            .split(|&b| b == b'\n')
            .filter(|l| l.starts_with(kind) && l.get(3) == Some(&b'\t'))
        {
            let mut attributes = BTreeMap::new();
            for field in line.split(|&b| b == b'\t').skip(1) {
                if field.len() < 3 || field[2] != b':' {
                    return Err(io::Error::other("invalid header attribute"));
                }
                attributes.insert(field[..2].to_vec(), field[3..].to_vec());
            }
            let id = attributes
                .remove(b"ID".as_slice())
                .ok_or_else(|| io::Error::other("header record without ID"))?;
            if !seen.insert(id.clone()) {
                return Err(io::Error::other("duplicate header ID"));
            }
            pending.push(Entry {
                file,
                id,
                attributes,
            });
        }
    }
    let pg = kind == b"@PG";
    let mut ready = Vec::new();
    pending.retain(|p| {
        if !pg || !p.attributes.contains_key(b"PP".as_slice()) {
            ready.push(p.clone());
            false
        } else {
            true
        }
    });
    let mut maps = vec![HashMap::new(); headers.len()];
    let mut taken = HashSet::new();
    let mut counter = 0;
    while !ready.is_empty() {
        let mut groups: RecordGroups = Vec::new();
        for p in &ready {
            let at = groups
                .iter()
                .position(|(id, _)| *id == p.id)
                .unwrap_or_else(|| {
                    groups.push((p.id.clone(), Vec::new()));
                    groups.len() - 1
                });
            let entries = &mut groups[at].1;
            if let Some((_, files)) = entries.iter_mut().find(|(attrs, _)| *attrs == p.attributes) {
                files.push(p.file);
            } else {
                entries.push((p.attributes.clone(), vec![p.file]));
            }
        }
        for (id, entries) in groups {
            for (_, files) in entries {
                let new = if taken.insert(id.clone()) {
                    counter += 1;
                    id.clone()
                } else {
                    loop {
                        let mut new = id.clone();
                        new.push(b'.');
                        new.extend_from_slice(base36(counter).as_bytes());
                        counter += 1;
                        if taken.insert(new.clone()) {
                            break new;
                        }
                    }
                };
                for file in files {
                    maps[file].insert(id.clone(), new.clone());
                }
            }
        }
        for p in &mut ready {
            if let Some(id) = maps[p.file].get(&p.id) {
                p.id = id.clone();
            }
        }
        for p in &mut pending {
            if let Some(id) = maps[p.file].get(&p.id) {
                p.id = id.clone();
            }
            if let Some(pp) = p.attributes.get_mut(b"PP".as_slice())
                && let Some(id) = maps[p.file].get(pp)
            {
                *pp = id.clone();
            }
        }
        let mut next = Vec::new();
        pending.retain(|p| {
            if ready
                .iter()
                .any(|r| r.file == p.file && p.attributes.get(b"PP".as_slice()) == Some(&r.id))
            {
                next.push(p.clone());
                false
            } else {
                true
            }
        });
        ready = next;
    }
    if !pending.is_empty() {
        return Err(io::Error::other("unresolved PG ancestry"));
    }
    Ok(maps)
}
