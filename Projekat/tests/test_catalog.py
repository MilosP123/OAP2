"""Provere postojećih notebook ćelija nad malim sintetičkim HDF5 katalogom."""
import json
import hashlib
import subprocess
import sys
from pathlib import Path

import h5py as hp
import numpy as np
import pandas as pd
import pytest

NOTEBOOK = Path(__file__).resolve().parents[1] / 'Untitled7.ipynb'
CELLS = json.loads(NOTEBOOK.read_text(encoding='utf-8'))['cells']


def run_cells(numbers, namespace):
    """Izvrši iste ćelije koje koristi notebook, bez kopiranja algoritama."""
    for number in numbers:
        if CELLS[number - 1]['cell_type'] != 'code':
            continue
        exec(compile(''.join(CELLS[number - 1]['source']),
                     f'notebook_cell_{number}', 'exec'), namespace)
    return namespace


def make_catalog(directory):
    """12 delova: leksikografski redosled mora pasti na broju 10."""
    directory.mkdir(parents=True)
    for chunk in range(12):
        with hp.File(directory / f'groupcat-99.{chunk}.hdf5', 'w') as data:
            header = data.create_group('Header')
            header.attrs.update({
                'NumFiles': 12, 'Ngroups_Total': 12, 'Nsubgroups_Total': 24,
                'Ngroups_ThisFile': 1, 'Nsubgroups_ThisFile': 2,
                'HubbleParam': 0.6774, 'Redshift': 0., 'Time': 1.,
                'BoxSize': 75000.,
                'SyntheticTest': 1,
            })
            group = data.create_group('Group')
            # Reprezentativne field/group/cluster mase.
            mass = [1e12, 1e13, 1e14][chunk % 3]
            group['Group_M_Crit200'] = [mass * 0.6774 / 1e10]
            group['GroupFirstSub'] = [2 * chunk]
            group['GroupNsubs'] = [2]
            sub = data.create_group('Subhalo')
            sub['SubhaloFlag'] = [1, 1]
            sub['SubhaloGrNr'] = [chunk, chunk]
            masses = np.zeros((2, 6))
            masses[:, 4] = [1., 2.]
            sub['SubhaloMassType'] = masses
            sub['SubhaloBHMass'] = [0.001 * (chunk + 1), 0.002]
            sub['SubhaloBHMdot'] = [0.002, 0.]
            sub['SubhaloSFR'] = [1., 0.]
            sub['SubhaloPos'] = np.full((2, 3), float(chunk))


@pytest.fixture
def working_catalog(tmp_path, monkeypatch):
    """Svaki test dobija svoje fajlove; nema pristupa originalnim podacima."""
    monkeypatch.chdir(tmp_path)
    make_catalog(tmp_path / 'data' / 'group99')
    return tmp_path


def pipeline():
    """Učitaj i izračunaj katalog istim redosledom kao u notebook-u."""
    return run_cells(range(1, 18), {})


def test_numeric_order_mapping_centrals_and_units(working_catalog):
    ns = pipeline()
    catalog = ns['clean']
    assert [int(Path(p).name.split('.')[-2]) for p in ns['files']] == list(range(12))
    np.testing.assert_array_equal(catalog['parent_group_id'], np.repeat(np.arange(12), 2))
    np.testing.assert_array_equal(catalog['is_central'], np.tile([True, False], 12))
    np.testing.assert_allclose(catalog['M200c_parent'], np.repeat([1e12, 1e13, 1e14] * 4, 2))
    assert catalog.iloc[0]['Mstar'] == pytest.approx(1e10 / .6774)
    assert catalog.iloc[0]['Mdot_BH'] == pytest.approx(.002 * 1e10 / .978e9)
    assert catalog.iloc[1]['lambda_edd'] == 0
    assert not catalog.iloc[1]['is_agn']
    assert len(catalog) == 24


def test_missing_chunk_fails(working_catalog):
    (working_catalog / 'data/group99/groupcat-99.4.hdf5').unlink()
    with pytest.raises(ValueError, match='Nedostaju delovi'):
        pipeline()


def test_missing_nonempty_field_fails(working_catalog):
    with hp.File(working_catalog / 'data/group99/groupcat-99.4.hdf5', 'a') as data:
        del data['Subhalo/SubhaloBHMdot']
    with pytest.raises(KeyError, match='SubhaloBHMdot'):
        pipeline()


def test_negative_parent_fails(working_catalog):
    with hp.File(working_catalog / 'data/group99/groupcat-99.4.hdf5', 'a') as data:
        data['Subhalo/SubhaloGrNr'][0] = -1
    with pytest.raises(ValueError, match='SubhaloGrNr'):
        pipeline()


def test_corrupted_central_mapping_fails(working_catalog):
    with hp.File(working_catalog / 'data/group99/groupcat-99.4.hdf5', 'a') as data:
        data['Group/GroupFirstSub'][0] = 0
    with pytest.raises(ValueError, match='GroupFirstSub'):
        pipeline()


def test_unknown_bh_values_and_zero_accretion(working_catalog):
    ns = pipeline()
    ns['MBH'][:6] = [0., 1e7, np.nan, np.inf, 1e7, 0.]
    ns['Mdot_BH'][:6] = [0., 0., .01, .01, -1., .01]
    ns['Mstar'][6] = np.inf
    run_cells([13, 16, 17], ns)
    status = ns['catalog']['is_agn']
    assert not status.iloc[0]       # no BH, no accretion
    assert not status.iloc[1]       # BH exists, zero accretion
    assert status.iloc[2:6].isna().all()
    assert 6 not in ns['clean'].index
    assert np.isnan(ns['catalog'].iloc[3]['logMBH'])
    summary = ns['summary']
    np.testing.assert_array_equal(summary['broj_nepoznatih'],
                                  summary['broj_galaksija'] - summary['broj_sa_AGN_statusom'])


def test_environment_boundaries():
    ns = {'np': np, 'M200c_parent': np.array([5e12 - 1, 5e12, 8e13 - 1, 8e13])}
    run_cells([5, 15], ns)
    assert ns['okruzenje'].tolist() == ['field', 'group', 'group', 'cluster-like']


def test_all_cells_and_export(working_catalog):
    import matplotlib
    matplotlib.use('Agg')
    ns = run_cells(range(1, 34), {})
    exported = pd.read_csv(working_catalog / 'outputs/catalog_base.csv.gz')
    assert len(exported) == len(ns['clean']) == 24
    assert (working_catalog / 'outputs/catalog_metadata.json').is_file()
    ns['plt'].close('all')


@pytest.mark.parametrize('field,values', [
    ('SubhaloPos', np.zeros((2, 2))),
    ('SubhaloMassType', np.zeros((2, 5))),
    ('SubhaloBHMass', np.zeros((2, 1))),
    ('SubhaloGrNr', np.array([4.1, 4.])),
])
def test_invalid_shapes_or_index_dtype(working_catalog, field, values):
    with hp.File(working_catalog / 'data/group99/groupcat-99.4.hdf5', 'a') as data:
        del data[f'Subhalo/{field}']
        data[f'Subhalo/{field}'] = values
    with pytest.raises(ValueError, match='oblik|celobrojno'):
        pipeline()


def test_satellite_cannot_be_mislabeled_central(working_catalog):
    with hp.File(working_catalog / 'data/group99/groupcat-99.4.hdf5', 'a') as data:
        data['Group/GroupFirstSub'][0] = 9  # Same halo, but its SECOND subhalo.
    with pytest.raises(ValueError, match='prvi globalni'):
        pipeline()


@pytest.mark.parametrize('key,value', [('HubbleParam', 0.), ('Time', .5),
                                      ('Ngroups_ThisFile', 1.5), ('BoxSize', np.inf)])
def test_bad_header(working_catalog, key, value):
    with hp.File(working_catalog / 'data/group99/groupcat-99.4.hdf5', 'a') as data:
        data['Header'].attrs[key] = value
    with pytest.raises(ValueError):
        pipeline()


def test_missing_header_field(working_catalog):
    with hp.File(working_catalog / 'data/group99/groupcat-99.4.hdf5', 'a') as data:
        del data['Header'].attrs['Time']
    with pytest.raises(KeyError, match='Header/Time'):
        pipeline()


def test_duplicate_chunk(working_catalog):
    import shutil
    directory = working_catalog / 'data/group99'
    shutil.copyfile(directory / 'groupcat-99.4.hdf5', directory / 'fof_subhalo_tab_099.4.hdf5')
    with pytest.raises(ValueError, match='duplikati'):
        pipeline()


def test_empty_chunk_may_omit_datasets(working_catalog):
    directory = working_catalog / 'data/group99'
    for path in directory.glob('*.hdf5'):
        with hp.File(path, 'a') as data:
            data['Header'].attrs['NumFiles'] = 13
            attrs = dict(data['Header'].attrs)
    with hp.File(directory / 'groupcat-99.12.hdf5', 'w') as data:
        attrs['Ngroups_ThisFile'] = 0
        attrs['Nsubgroups_ThisFile'] = 0
        data.create_group('Header').attrs.update(attrs)
    assert len(pipeline()['clean']) == 24


@pytest.mark.parametrize('unsigned', [False, True])
def test_group_without_subhalos(working_catalog, unsigned):
    directory = working_catalog / 'data/group99'
    for path in directory.glob('*.hdf5'):
        with hp.File(path, 'a') as data:
            data['Header'].attrs['Ngroups_Total'] = 13
    with hp.File(directory / 'groupcat-99.11.hdf5', 'a') as data:
        data['Header'].attrs['Ngroups_ThisFile'] = 2
        for field, values in {
            'Group_M_Crit200': [1e4, 1.], 'GroupNsubs': [2, 0],
            'GroupFirstSub': np.array([22, 2**32 - 1], dtype=np.uint32) if unsigned
            else np.array([22, -1], dtype=np.int32),
        }.items():
            del data[f'Group/{field}']
            data[f'Group/{field}'] = values
    ns = pipeline()
    assert ns['is_central'].sum() == 12
    assert ns['GroupFirstSub'][-1] == -1


def test_eddington_normalization_and_agn_threshold(working_catalog):
    ns = pipeline()
    # Independent reference: 0.1 Msun/year at 1e8 Msun gives lambda about 0.04494.
    ns['MBH'][:] = 1e8
    ns['Mdot_BH'][:] = .1
    run_cells([13, 16, 17], ns)
    assert ns['lambda_edd'][0] == pytest.approx(.04494, rel=2e-4)
    assert ns['clean']['is_agn'].all()
    ns['lambda_edd'][:3] = [1e-3 - 1e-10, 1e-3, 1e-3 + 1e-10]
    run_cells([16], ns)
    assert ns['catalog']['is_agn'].iloc[:3].tolist() == [False, True, True]


def test_unknown_status_denominator_and_bounds(working_catalog):
    ns = pipeline()
    frame = ns['clean'].iloc[:3].copy()
    frame['is_agn'] = pd.array([True, False, pd.NA], dtype='boolean')
    frame['environment'] = 'field'
    result = ns['pregled_uzorka'](frame, 'environment').iloc[0]
    assert result['broj_galaksija'] == 3
    assert result['broj_sa_AGN_statusom'] == 2
    assert result['AGN_frakcija'] == .5
    assert result['AGN_frakcija_min'] == pytest.approx(1/3)
    assert result['AGN_frakcija_max'] == pytest.approx(2/3)


@pytest.mark.parametrize('empty_base', [False, True])
def test_no_positive_accretion_and_empty_sample(working_catalog, empty_base):
    for path in (working_catalog / 'data/group99').glob('*.hdf5'):
        with hp.File(path, 'a') as data:
            data['Subhalo/SubhaloBHMdot'][:] = 0.
            if empty_base:
                data['Subhalo/SubhaloFlag'][:] = 0
    ns = run_cells(range(1, 34), {})
    assert len(ns['clean']) == (0 if empty_base else 24)
    assert ns['clean']['is_agn'].sum() == 0
    assert len(ns['plt'].get_fignums()) == 0


def test_halo_bootstrap_and_mass_bin_edges(working_catalog):
    ns = pipeline()
    ns['clean']['logMstar'] = 10.  # Exact bin boundary, including the maximum.
    run_cells([27], ns)
    assert ns['mass_fractions']['N'].sum() == len(ns['clean'])
    frame = ns['clean'].iloc[:3].copy()
    frame['parent_group_id'] = [0, 0, 1]
    frame['is_agn'] = pd.array([True, True, False], dtype='boolean')
    a = ns['halo_interval'](frame, np.random.default_rng(42), 1000)
    b = ns['halo_interval'](frame, np.random.default_rng(42), 1000)
    assert a == b
    assert a == (0., 1., 2)
    one_halo = ns['halo_interval'](frame.iloc[:2], np.random.default_rng(42))
    assert np.isnan(one_halo[0]) and one_halo[2] == 1


def test_input_bytes_unchanged_after_full_execution(working_catalog):
    paths = sorted((working_catalog / 'data/group99').glob('*.hdf5'))
    before = [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]
    ns = run_cells(range(1, 34), {})
    after = [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]
    assert before == after
    metadata = json.loads((working_catalog / 'outputs/catalog_metadata.json').read_text())
    assert metadata['input_kind'] == 'synthetic_test'
    exported = working_catalog / 'outputs/catalog_base.csv.gz'
    assert metadata['catalog_sha256'] == hashlib.sha256(exported.read_bytes()).hexdigest()
    assert not ns['plt'].get_fignums()


def test_fresh_python_process(working_catalog):
    code = '''import json, sys
from pathlib import Path
nb = json.loads(Path(sys.argv[1]).read_text())
namespace = {}
for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        exec(compile(''.join(cell['source']), '<notebook>', 'exec'), namespace)
assert len(namespace['clean']) == 24
'''
    result = subprocess.run([sys.executable, '-c', code, str(NOTEBOOK)],
                            cwd=working_catalog, capture_output=True, text=True,
                            timeout=60, check=False)
    assert result.returncode == 0, result.stderr


def test_missing_data_fails_clearly(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError, match='Nema delova kataloga'):
        pipeline()


def test_finite_halo_selection_and_unknown_environment(working_catalog):
    ns = pipeline()
    ns['M200c_parent'][:4] = [np.nan, np.inf, 0., -1.]
    run_cells([15, 16], ns)
    assert len(ns['clean']) == 20
    assert (ns['catalog']['environment'].iloc[:4] == 'unknown').all()


@pytest.mark.parametrize('value', [np.nan, np.inf, -np.inf, -1.])
def test_invalid_accretion_is_unknown(working_catalog, value):
    ns = pipeline()
    ns['Mdot_BH'][0] = value
    run_cells([13, 16], ns)
    assert pd.isna(ns['catalog'].iloc[0]['is_agn'])
    assert pd.isna(ns['catalog'].iloc[0]['lambda_edd'])


def test_positive_out_of_range_parent_rejected(working_catalog):
    with hp.File(working_catalog / 'data/group99/groupcat-99.4.hdf5', 'a') as data:
        data['Subhalo/SubhaloGrNr'][0] = 12
    with pytest.raises(ValueError, match='SubhaloGrNr'):
        pipeline()
