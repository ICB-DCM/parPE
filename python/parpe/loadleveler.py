from jinja2 import Template, Environment, PackageLoader

env = Environment(loader=PackageLoader('parpe', 'templates'), autoescape=True)


def create_job_file():
    num_starts = 3
    ll = {
        'job_name': 'testname',
        'steps': [],
        'input_data_file': 'input.h5',
    }

    ll['steps'].append({
        'step_name': 'preprocess',
        'class': 'micro',
        'node': 1,
        'wall_clock_limit': '02:00:00',
        'energy_policy_tag':
            'simulation2_energy_tag_ipopt_hierarchical_test',
        'body': 'snakemake preprocess'
    })

    for start_idx in range(num_starts):
        ll['steps'].append({
            'step_name': f'optimize_ms_{start_idx}',
            'class': 'micro',
            'node': 1,
            'wall_clock_limit': '02:00:00',
            'energy_policy_tag':
                'simulation2_energy_tag_ipopt_hierarchical_test',
            'dependency': '(preprocess == 0)',
            'body': 'snakemake optimize' # TODO: start idx
        })

    ll['steps'].append({
        'step_name': 'postprocess',
        'class': 'micro',
        'node': 1,
        'wall_clock_limit': '02:00:00',
        'energy_policy_tag':
            'simulation2_energy_tag_ipopt_hierarchical_test',
        'dependency': ' && '.join([f'(optimize_ms_{start_idx} >= -999)'
                                   for start_idx in range(num_starts)]),
        'body': 'snakemake postprocess'
    })

    template = env.get_template('loadleveler_petab.tpl')
    res = template.render(ll=ll)
    print(res)


if __name__ == '__main__':
    create_job_file()
