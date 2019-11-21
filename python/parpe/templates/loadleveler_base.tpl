#!/bin/bash
# LoadLeveler job file for running an optimization using parPE and running
# and running simulations at the optimal point.
#
# 2017-19 Daniel Weindl <daniel.weindl@helmholtz-muenchen.de>

# General settings
# @ job_name = {{ ll.get('job_name', 'Unnamed job') }}
{%- if ll.comment %}# @ comment = {{ll.comment}}{% endif %}
# @ output = {% if ll.output %}{{ ll.output }}{% else %}job-$(user)-$(jobid)-$(job_name)-$(stepid)-$(step_name).out{% endif %}
# @ error = {% if ll.error %}{{ ll.error }}{% else %}job-$(user)-$(jobid)-$(job_name)-$(stepid)-$(step_name).err{% endif %}
# @ island_count = 1
# @ minimize_time_to_solution = yes
# @ adjust_wall_clock_limit = yes
# @ notification=always
{%- if ll.notify_user %}# @ notify_user = {{ll.notify_user}}{% endif %}
# @ network.MPI = sn_all,not_shared,us

{% block body %}
{% endblock %}
