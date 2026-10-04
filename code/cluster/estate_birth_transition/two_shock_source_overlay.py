"""Exact-hash read-only overlay for legacy absolute source paths in the frozen observer."""
from __future__ import annotations
import builtins
import hashlib
import importlib.machinery
import io
import json
import os
import shutil
from pathlib import Path

EXPECTED_SOURCE_MANIFEST_SHA256 = '07d84336a3112b251afe505908113d9c00585b91c34bd50f0dee108435db496d'


def _sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def prepare_overlay(repo_root, frozen_root, v9_inventory):
    """Snapshot exact manifest bytes, recovering only through authenticated v9 inventory entries."""
    repo=Path(repo_root).resolve(); frozen=Path(frozen_root).resolve()
    manifest=repo/'output/model/overnight_calibration_20260928/contract_v1/source_manifest.json'
    if _sha(manifest)!=EXPECTED_SOURCE_MANIFEST_SHA256: raise RuntimeError('Nested source manifest pin mismatch')
    contract=json.loads(manifest.read_text()); files=contract['files']
    inv=json.loads(Path(v9_inventory).read_text())['files']
    candidates={}
    for rel,digest in inv.items(): candidates.setdefault(digest,[]).append(rel)
    base=frozen/'source_overlay'; snapshot=base/'files'; mapping={}; recovered=[]
    for rel in files: mapping[str((repo/rel).resolve())]=str((snapshot/rel).resolve())
    if base.exists():
        existing=base/'overlay.json'
        if not existing.is_file(): raise RuntimeError('Refusing to replace incomplete source overlay directory')
        prior=json.loads(existing.read_text())
        if prior.get('schema')!='exact_read_only_source_overlay_v1' or prior.get('manifest_sha256')!=EXPECTED_SOURCE_MANIFEST_SHA256 or prior.get('original_root')!=str(repo) or prior.get('mapping')!=mapping:
            raise RuntimeError('Existing source overlay identity differs; refusing replacement')
        for rel,digest in files.items():
            if _sha(snapshot/rel)!=digest: raise RuntimeError(f'Existing overlay snapshot hash mismatch: {rel}')
            (snapshot/rel).chmod(0o444)
        return prior
    for rel,digest in files.items():
        original=repo/rel
        source=None; source_rel=None
        if original.is_file() and _sha(original)==digest:
            source=original; source_rel='workspace-exact'
        else:
            exact=[candidate for candidate in candidates.get(digest,[]) if (frozen/'source'/candidate).is_file() and _sha(frozen/'source'/candidate)==digest]
            if not exact: raise RuntimeError(f'No exact authenticated bytes for nested source pin: {rel} sha256={digest}')
            source=frozen/'source'/exact[0]; source_rel='deployment_v9:'+exact[0]
            recovered.append({'relative_path':rel,'sha256':digest,'source':source_rel})
        dest=snapshot/rel; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(source,dest)
        if _sha(dest)!=digest: raise RuntimeError(f'Overlay copy hash mismatch: {rel}')
        dest.chmod(0o444)
    metadata={'schema':'exact_read_only_source_overlay_v1','manifest_path':str(manifest),
        'manifest_sha256':EXPECTED_SOURCE_MANIFEST_SHA256,'manifest_file_count':len(files),
        'original_root':str(repo),'snapshot_root':str(snapshot),'mapping':mapping,'recovered_from_v9':recovered,
        'all_source_hashes_verified':True}
    mp=base/'overlay.json'; mp.write_text(json.dumps(metadata,indent=2,sort_keys=True)+'\n')
    return metadata


def install_overlay(overlay_json):
    """Map exact legacy paths to read-only snapshots, retaining original import filenames."""
    import sys
    marker='_two_shock_exact_source_overlay'
    if getattr(sys,marker,None): raise RuntimeError('Source overlay already installed')
    meta=json.loads(Path(overlay_json).read_text())
    if meta.get('schema')!='exact_read_only_source_overlay_v1' or not meta.get('all_source_hashes_verified'):
        raise RuntimeError('Authenticated source overlay receipt required')
    if meta.get('manifest_sha256')!=EXPECTED_SOURCE_MANIFEST_SHA256:
        raise RuntimeError('Source overlay manifest is not the pinned contract')
    mapping={str(Path(a).resolve()):Path(b).resolve() for a,b in meta['mapping'].items()}
    original_root=Path(meta.get('original_root','')).resolve()
    snapshot_root=Path(meta['snapshot_root']).resolve()
    source_manifest=Path(meta['manifest_path'])
    if _sha(source_manifest)!=meta['manifest_sha256']: raise RuntimeError('Source manifest changed before install')
    expected=json.loads(source_manifest.read_text())['files']
    expected_mapping={str((original_root/rel).resolve()): (snapshot_root/rel).resolve() for rel in expected}
    if mapping!=expected_mapping: raise RuntimeError('Overlay mapping differs from exact manifest roots/paths')
    if len(mapping)!=meta.get('manifest_file_count') or len(mapping)!=len(expected): raise RuntimeError('Overlay file count differs from manifest')
    for rel,digest in expected.items():
        target=expected_mapping[str((original_root/rel).resolve())]
        if _sha(target)!=digest: raise RuntimeError('Overlay snapshot differs from source contract: '+rel)
        if target.stat().st_mode & 0o222: raise RuntimeError('Overlay snapshot has write permission: '+rel)
    mapped_paths=set(mapping)|{str(p) for p in mapping.values()}
    def resolve(file):
        try: return mapping.get(str(Path(file).resolve()),file)
        except (TypeError,ValueError,OSError): return file
    def mode_write(mode): return any(flag in mode for flag in ('w','a','x','+'))
    old_builtin=builtins.open; old_io=io.open; old_path=Path.open
    old_get_data=importlib.machinery.SourceFileLoader.get_data
    old_get_code=importlib.machinery.SourceFileLoader.get_code
    def safe_builtin(file,mode='r',*args,**kwargs):
        if isinstance(file,(str,Path)) and str(Path(file).resolve()) in mapped_paths and mode_write(mode): raise PermissionError('Authenticated source overlay is read-only')
        return old_builtin(resolve(file),mode,*args,**kwargs)
    def safe_io(file,mode='r',*args,**kwargs):
        if isinstance(file,(str,Path)) and str(Path(file).resolve()) in mapped_paths and mode_write(mode): raise PermissionError('Authenticated source overlay is read-only')
        return old_io(resolve(file),mode,*args,**kwargs)
    def safe_get_data(loader,path): return old_get_data(loader,str(resolve(path)))
    def safe_get_code(loader,fullname):
        original=str(Path(loader.path).resolve())
        if original in mapping and original.endswith('.py'):
            return loader.source_to_code(mapping[original].read_bytes(),loader.path)
        return old_get_code(loader,fullname)
    old_os_open=os.open
    write_flags=os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND
    def safe_os_open(path,flags,*args,**kwargs):
        if isinstance(path,(str,Path)) and str(Path(path).resolve()) in mapped_paths and flags & write_flags:
            raise PermissionError('Authenticated source overlay is read-only')
        return old_os_open(resolve(path),flags,*args,**kwargs)
    def safe_path(self,mode='r',*args,**kwargs):
        if str(Path(self).resolve()) in mapped_paths and mode_write(mode): raise PermissionError('Authenticated source overlay is read-only')
        return old_path(self if mode_write(mode) else resolve(self),mode,*args,**kwargs)
    builtins.open=safe_builtin; io.open=safe_io; Path.open=safe_path; os.open=safe_os_open
    importlib.machinery.SourceFileLoader.get_data=safe_get_data
    importlib.machinery.SourceFileLoader.get_code=safe_get_code
    receipt={'manifest_sha256':meta['manifest_sha256'],'mapped_files':len(mapping),'snapshot_root':meta['snapshot_root'],
        'write_denial':'builtins.open, io.open, Path.open reject write modes for mapped source paths',
        'loader_mapping':'SourceFileLoader.get_code compiles exact snapshot bytes, bypassing cached bytecode; original loader path and __file__ are retained'}
    setattr(sys,marker,receipt)
    return receipt
