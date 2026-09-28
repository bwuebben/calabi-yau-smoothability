#!/usr/bin/env python3
"""Check bundle file integrity. This does not replace the mathematical checkers."""
import hashlib,json
from pathlib import Path


def verify(root):
    manifest=json.loads((root/'MANIFEST.json').read_text())
    assert manifest['schema_version']==1 and manifest['status']=='complete'
    for rel,entry in manifest['files'].items():
        path=root/rel
        if rel=='README.md' and (root/'COMPUTATIONS.md').exists():path=root/'COMPUTATIONS.md'
        assert path.is_file(),f'Missing bundle file: {rel}'
        data=path.read_bytes()
        assert len(data)==entry['bytes'] and hashlib.sha256(data).hexdigest()==entry['sha256'],f'Changed bundle file: {rel}'
    for key in ['x9','x19','x20','x88','x154']:assert f'examples/{key}.json' in manifest['files']
    assert len([p for p in manifest['files'] if p.startswith('unimodular/certificates/')])==30
    print(f"PASS: {len(manifest['files'])} bundle files match the manifest.")
    return manifest

if __name__=='__main__':verify(Path(__file__).resolve().parent)
