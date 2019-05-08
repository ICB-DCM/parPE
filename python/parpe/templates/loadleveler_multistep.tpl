{% extends "loadleveler_base.tpl" %}
{%- block body %}
    {%- for step in ll.steps %}
# Step {{ loop.index }}
{% include "loadleveler_step_header.tpl" %}
    {% endfor %}

{%- block before_any_step %}{{ ll.get('before_any_step', '') }}{% endblock %}

case ${LOADL_STEP_NAME} in
    {%- for step in ll.steps %}
    {{ step.step_name }})
        {{ step.body }}
    ;;
    {% endfor %}
esac

{%- block after_any_step %}{{ ll.get('after_any_step', '') }}{% endblock %}

{% endblock %}
