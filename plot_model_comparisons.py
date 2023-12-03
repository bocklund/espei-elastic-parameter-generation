import pathlib
import matplotlib.pyplot as plt
import tinydb
import pycalphad
import numpy as np
from pycalphad import Database, calculate, Model
from espei.datasets import load_datasets, recursive_glob
from espei.core_utils import get_prop_data, filter_configurations, filter_temperatures, symmetry_filter
from espei.sublattice_tools import canonicalize

class ElasticModel(Model):
    def build_phase(self, dbe):
        super().build_phase(dbe)
        phase = dbe.phases[self.phase_name]
        param_search = dbe.search
        for prop in ['C11', 'C12', 'C44']:
            prop_param_query = (
                (tinydb.where('phase_name') == phase.name) & \
                (tinydb.where('parameter_type') == prop) & \
                (tinydb.where('constituent_array').test(self._array_validity))
                )
            prop_val = self.redlich_kister_sum(phase, param_search, prop_param_query).subs(dbe.symbols)
            setattr(self, prop, prop_val)


def compare_data_to_model(dbf, comps, phase_name, configuration, output_property, datasets, x_axis_component, ax):
    solver_qry = (tinydb.where('solver').test(symmetry_filter, configuration, None))
    desired_data = get_prop_data(comps, phase_name, [output_property], datasets, additional_query=solver_qry)
    desired_data = filter_configurations(desired_data, configuration, None)
    desired_data = filter_temperatures(desired_data)
    # we know that there's only one dataset and they are directly the absolute values with no reference state to worry about

    ds = desired_data[0]
    exp_x_data = [occup[0][config[0].index(x_axis_component)] for config, occup in zip(ds["solver"]["sublattice_configurations"], ds["solver"]["sublattice_occupancies"])]  # assumption about the sublattice model here
    exp_y_data = ds["values"]

    x_comp_config_idx = canonicalize(configuration, None)[0].index(x_axis_component)
    num_points = 1001
    calc_x_data = np.linspace(0, 1, num_points)
    # assumption about the configuration being ((A, B,), VA)
    points = np.empty((num_points, 3))
    points[:, x_comp_config_idx] = calc_x_data
    points[:, 1 - x_comp_config_idx] = 1 - calc_x_data
    points[:, 2] = np.ones_like(calc_x_data)
    calc_res = calculate(dbf, comps, phase_name, N=1, P=101325, T=298.15, output=output_property, points=points, model=ElasticModel)
    calc_y_data = calc_res[output_property].squeeze()

    ax.plot(calc_x_data, calc_y_data)
    ax.scatter(exp_x_data, exp_y_data)
    ax.set_xlabel(f"X({x_axis_component})")
    ax.set_xlim(0, 1)
    ax.set_ylabel(output_property)
    ax.set_title(f"{output_property}: {configuration}")


# PyCalphad's TDB reader doesn't know about our custom parameters, so we need to add them ad valid options
pycalphad.io.tdb_keywords.TDB_PARAM_TYPES.extend(['C11', 'C12', 'C44'])

# now we can load our databases
dbf = Database("Ti-elastic.tdb")
# and datasets
datasets = load_datasets(recursive_glob("elastic-datasets", '*.json'))

outdir = pathlib.Path("figures")
outdir.mkdir(parents=True, exist_ok=True)

output_property = "C11"
fig, (ax_TiMo, ax_TiZr) = plt.subplots(ncols=2, figsize=(10, 4))
compare_data_to_model(dbf, ["MO", "TI", "VA"], "BCC_A2", (('MO', 'TI'), "VA"), output_property, datasets, "MO", ax_TiMo)
compare_data_to_model(dbf, ["TI", "ZR", "VA"], "BCC_A2", (('TI', 'ZR'), "VA"), output_property, datasets, "ZR", ax_TiZr)
ax_TiMo.set_ylim(0, 500)
ax_TiZr.set_ylim(0, 200)
filename = outdir / f"{output_property}.png"
fig.savefig(filename, dpi=300, bbox_inches="tight")

output_property = "C12"
fig, (ax_TiMo, ax_TiZr) = plt.subplots(ncols=2, figsize=(10, 4))
compare_data_to_model(dbf, ["MO", "TI", "VA"], "BCC_A2", (('MO', 'TI'), "VA"), output_property, datasets, "MO", ax_TiMo)
compare_data_to_model(dbf, ["TI", "ZR", "VA"], "BCC_A2", (('TI', 'ZR'), "VA"), output_property, datasets, "ZR", ax_TiZr)
ax_TiMo.set_ylim(0, 200)
ax_TiZr.set_ylim(0, 140)
filename = outdir / f"{output_property}.png"
fig.savefig(filename, dpi=300, bbox_inches="tight")

output_property = "C44"
fig, (ax_TiMo, ax_TiZr) = plt.subplots(ncols=2, figsize=(10, 4))
compare_data_to_model(dbf, ["MO", "TI", "VA"], "BCC_A2", (('MO', 'TI'), "VA"), output_property, datasets, "MO", ax_TiMo)
compare_data_to_model(dbf, ["TI", "ZR", "VA"], "BCC_A2", (('TI', 'ZR'), "VA"), output_property, datasets, "ZR", ax_TiZr)
ax_TiMo.set_ylim(0, 120)
ax_TiZr.set_ylim(0, 70)
filename = outdir / f"{output_property}.png"
fig.savefig(filename, dpi=300, bbox_inches="tight")
