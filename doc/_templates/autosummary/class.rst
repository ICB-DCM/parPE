{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :show-inheritance:
   :inherited-members:

   {% block methods %}
   {% if '__init__' in methods %}
   .. automethod:: __init__
      :noindex:
   {% endif %}
   {% if methods %}{% if functions %}
.. rubric:: Functions

{% for function in functions %}
.. autofunction:: {{ function }}
{% endfor %}

{% endif %}

   .. rubric:: Methods Summary
   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes
   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods_full %}
   {% if methods %}
   .. rubric:: Methods
   {% for item in methods %}
   .. automethod:: {{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}


