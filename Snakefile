configfile: 'config/config.yaml'

from workflow.config_adapter import merge_cwl_job_yaml_into_config


if config.get('cwl_input_yaml'):
    config = merge_cwl_job_yaml_into_config(config, config['cwl_input_yaml'])

include: 'workflow/rules/se_foundation.smk'
