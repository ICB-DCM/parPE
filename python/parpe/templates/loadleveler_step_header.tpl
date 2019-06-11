# @ step_name = {{ step.get('step_name', 'unnamed_step') }}
{%- if step.dependency %}# @ dependency = {{step.dependency}}{% endif %}
# @ class = {{ step.get('class', 'test') }}
# @ node = {{ step.get('node', 1) }}
# @ tasks_per_node = {{ step.get('tasks_per_node', 28) }}
# @ wall_clock_limit = {{ step.get('wall_clock_limit', '00:30:00') }}
# @ job_type = parallel
# @ energy_policy_tag = {{ step.get('energy_policy_tag', '') }}
# @ queue
