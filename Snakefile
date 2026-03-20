configfile: 'config/config.yaml'

from workflow.config_adapter import merge_cwl_job_yaml_into_config


if config.get('cwl_input_yaml'):
    config = merge_cwl_job_yaml_into_config(config, config['cwl_input_yaml'])

pipeline_profile = config.get('pipeline_profile', 'mini')

if pipeline_profile == 'mini':
    include: 'workflow/rules/se_foundation.smk'
elif pipeline_profile == 'dropin':
    include: 'workflow/rules/dropin_repelement.smk'
else:
    raise ValueError(f'Unknown pipeline_profile={pipeline_profile!r}; expected mini or dropin')
