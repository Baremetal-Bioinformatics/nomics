"""Small local script checks; run with python3 -I scripts/goldens/tests/somatic_scripts_test.py."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[3]


def run(*args, **kwargs):
    return subprocess.run(args, text=True, capture_output=True, **kwargs)


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    pipeline = root / 'run/results/pipeline_info'
    pipeline.mkdir(parents=True)
    source = root / 'inputs'
    source.mkdir()
    (source / 'a.bed').write_text('chr1\t0\t10\n')
    (source / 'reads.cram').write_text('mock cram')
    rows = ['task_id\thash\tnative_id\tname\tstatus\n']
    for i, status in enumerate(['COMPLETED', 'FAILED', 'COMPLETED'], 1):
        work = root / f'run/work/aa/{i}abc-full'
        work.mkdir(parents=True)
        (work / 'a.bed').symlink_to(source / 'a.bed')
        (work / 'reads.cram').symlink_to(source / 'reads.cram')
        (work / '.command.sh').write_text('gatk -I reads.cram -L a.bed\n')
        (work / 'output.vcf').write_text(str(i))
        rows.append(f'{i}\taa/{i}abc\t0\tA:BAM_VARIANT_CALLING_SOMATIC_MUTECT2:MUTECT2_PAIRED (tumor_vs_normal)\t{status}\n')
    (pipeline / 'execution_trace_test.txt').write_text(''.join(rows))
    harvest = REPO / 'scripts/goldens/somatic_harvest.sh'
    result = run('bash', str(harvest), str(root / 'run'))
    assert result.returncode == 0, result.stderr
    golden = root / 'run/golden'
    assert (golden / 'MUTECT2_PAIRED-tumor_vs_normal.2/output.vcf').read_text() == '3'
    assert not (golden / 'MUTECT2_PAIRED-tumor_vs_normal/reads.cram').exists()
    assert 'reads.cram -> ' in (golden / 'MUTECT2_PAIRED-tumor_vs_normal/INPUTS').read_text()
    for line in (golden / 'SHA256SUMS').read_text().splitlines():
        digest, path = line.split('  ', 1)
        assert hashlib.sha256((golden / path).read_bytes()).hexdigest() == digest
    sums = (golden / 'SHA256SUMS').read_bytes()
    assert run('bash', str(harvest), str(root / 'run')).returncode != 0
    assert run('bash', str(harvest), str(root / 'run'), '--force').returncode == 0
    assert sums == (golden / 'SHA256SUMS').read_bytes()
    print('PASS mock harvest: trace order, failed-task exclusion, CRAM provenance, BED, hashes, overwrite refusal, deterministic rebuild')

    # Prepare tiny fake inputs and check coordinate/order validation and BAM cache reuse.
    data = root / 'prepared'
    art = root / 'art'
    data.mkdir(); art.mkdir()
    dictionary = source / 'ref.dict'
    dictionary.write_text('@SQ\tSN:chr2\tLN:100\n@SQ\tSN:chr1\tLN:100\n')
    intervals = source / 'targets.interval_list'
    intervals.write_text('@HD\tVN:1.6\nchr2\t1\t10\t+\ta\nchr1\t11\t20\t+\tb\n')
    env = dict(os.environ, REPO=str(REPO), DATA=str(data), ART=str(art), DICT=str(dictionary), INTERVALS=str(intervals), BUNDLE=str(source), PATIENT='HCC1143', SEX='XX', TUMOUR_SAMPLE='tumor', NORMAL_SAMPLE='normal')
    for role in ('TUMOUR', 'NORMAL'):
        for kind in ('BAM', 'BAI'):
            p = source / f'{role}.{kind.lower()}'
            p.write_text('test')
            env[role + '_' + kind] = str(p)
    for name in ('af-only-gnomad.hg38.vcf.gz', 'af-only-gnomad.hg38.vcf.gz.tbi', '1000g_pon.hg38.vcf.gz', '1000g_pon.hg38.vcf.gz.tbi'):
        (source / name).write_text('test')
    def prepare():
        return run('python3', '-I', str(REPO / 'scripts/goldens/somatic_bam_prepare.py'), env=env)
    result = prepare()
    assert result.returncode == 0, result.stderr
    assert (data / 'intervals.bed').read_text() == 'chr2\t0\t10\nchr1\t10\t20\n'
    params = (data / 'params.yaml').read_text()
    assert 'step: "variant_calling"' in params and 'tools: "mutect2"' in params
    assert 'save_output_as_bam: true' in params and 'publish_dir_mode: link' in params
    assert not any(params.startswith(k + ':') or '\n' + k + ':' in params for k in ('bwa', 'aligner', 'known_indels', 'skip_tools'))
    sheet = (data / 'samplesheet.csv').read_text()
    assert sheet.startswith('patient,sex,status,sample,bam,bai\nHCC1143,XX,0,normal,')
    assert '\nHCC1143,XX,1,tumor,' in sheet
    cache = Path(env['TUMOUR_BAM'] + '.md5.json')
    mtime = cache.stat().st_mtime_ns
    assert prepare().returncode == 0 and cache.stat().st_mtime_ns == mtime
    for bad in ('chr1\t1\t10\nchr2\t1\t10\n', 'chr2\t1\t10\nchr2\t10\t20\n', 'chrX\t1\t10\n', 'chr2\t0\t10\n'):
        intervals.write_text(bad)
        assert prepare().returncode != 0
    print('PASS preparation: Picard conversion, dictionary order, overlap/bounds/unknown-contig rejection, checksum cache reuse')

    # Capture all nomics argv from old/new gates, using local S1 goldens only.
    fixture = Path.home() / 'omics/runs/somatic/S1-r1/golden'
    if fixture.is_dir():
        gate_run = root / 'gate-run'
        gate_run.mkdir()
        (gate_run / 'golden').symlink_to(fixture)
        stub = root / 'nomics'
        stub.write_text('#!/usr/bin/env python3\nimport json, os, sys\nwith open(os.environ["ARGV_LOG"], "a") as f: f.write(json.dumps(sys.argv[1:])+"\\n")\nsys.exit(1)\n')
        stub.chmod(0o755)
        logs = {}
        for version in ('before', 'after'):
            for name in ('somatic_boundary1_gate.sh', 'somatic_chain_gate.sh'):
                if version == 'before':
                    text = run('git', 'show', f'ea9403f:{"scripts/" + name}', cwd=REPO).stdout
                else:
                    text = (REPO / 'scripts' / name).read_text()
                text = text.replace('REPO=$(cd "$(dirname "$0")/.." && pwd)', f'REPO="{REPO}"')
                text = text.replace('/usr/bin/time -v ', '')  # GNU time is unavailable on this Mac.
                # Test-only compatibility for the Linux gate's mapfile -t call on Bash 3.
                text = text.replace('set -uo pipefail', r'''set -uo pipefail
mapfile() { local line; eval "$2=()"; while IFS= read -r line; do eval "$2+=(\"\$line\")"; done; }
''')
                text = text.replace('$REPO/scripts/somatic_boundary1_gate.sh', str(root / f'{version}-somatic_boundary1_gate.sh'))
                (root / f'{version}-{name}').write_text(text)
            for mode in ('steps', 'one', 'boundary'):
                log = root / f'{version}-{mode}.jsonl'
                out = root / 'gate-out'
                shutil.rmtree(out, ignore_errors=True)
                env2 = dict(os.environ, NOMICS=str(stub), ARGV_LOG=str(log), CHAIN_ONE='1' if mode == 'one' else '0')
                for key in ('CHAIN_REUSE', 'TUMOUR_SAMPLE', 'NORMAL_SAMPLE', 'NORMAL_SM', 'TUMOUR_BAM', 'NORMAL_BAM', 'TUMOUR_DIR', 'NORMAL_DIR'):
                    env2.pop(key, None)
                script = 'somatic_boundary1_gate.sh' if mode == 'boundary' else 'somatic_chain_gate.sh'
                result = run('bash', str(root / f'{version}-{script}'), str(gate_run), str(out), env=env2)
                assert log.exists(), result.stderr
                logs[version, mode] = log.read_text()
        for mode in ('steps', 'one', 'boundary'):
            assert logs['before', mode] == logs['after', mode], mode
            print(f'PASS default argv diff ({mode}): {len(logs["after", mode].splitlines())} invocations, zero differences')
        custom = root / 'custom-run'
        custom_golden = custom / 'golden'
        custom_golden.mkdir(parents=True)
        (custom_golden / 'intervals').symlink_to(fixture / 'intervals')
        for directory in fixture.iterdir():
            if directory.is_dir() and (directory / '.command.sh').exists():
                label = directory.name.replace('HCC1395BL_WES', 'normal').replace('HCC1395_WES', 'tumor')
                dest = custom_golden / label
                dest.mkdir()
                command = (directory / '.command.sh').read_text().replace('HCC1395BL_WES', 'normal').replace('HCC1395_WES', 'tumor')
                (dest / '.command.sh').write_text(command)
        log = root / 'override.jsonl'
        env2.update(ARGV_LOG=str(log), CHAIN_ONE='0', TUMOUR_SAMPLE='tumor', NORMAL_SAMPLE='normal', NORMAL_SM='HCC1143_normal', TUMOUR_BAM='/fixture/tumor.bam', NORMAL_BAM='/fixture/normal.bam')
        run('bash', str(root / 'after-somatic_chain_gate.sh'), str(custom), str(root / 'custom-out'), env=env2)
        calls = [json.loads(line) for line in log.read_text().splitlines()]
        mutect = next(c for c in calls if c[:2] == ['somatic', 'mutect2'])
        assert '/fixture/normal.bam' in mutect and '/fixture/tumor.bam' in mutect
        assert mutect[mutect.index('--normal-sample') + 1] == 'HCC1143_normal'
        pileups = [c for c in calls if c[:2] == ['somatic', 'get-pileup-summaries']]
        assert len(pileups) == 4 and all('/fixture/' in c[c.index('--input') + 1] for c in pileups)
        env2.update(ARGV_LOG=str(root / 'override-one.jsonl'), CHAIN_ONE='1')
        run('bash', str(root / 'after-somatic_chain_gate.sh'), str(custom), str(root / 'custom-one'), env=env2)
        chain = json.loads((root / 'override-one.jsonl').read_text())
        assert chain[chain.index('--normal') + 1] == '/fixture/normal.bam'
        assert chain[chain.index('--normal-label') + 1] == 'normal'
        assert chain[chain.index('--tumour-label') + 1] == 'tumor'
        print('PASS custom BAM paths, labels and normal SM in both gate modes and nested boundary gate')
    else:
        raise RuntimeError('S1 goldens required for default argv differential check')
