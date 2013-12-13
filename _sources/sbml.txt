.. _section-sbml:

Inchman SBML Annotations
========================
The current version of SBML does not yet support spatial information. For this reason, Inchman uses an annotated SBML model to convey any spatial and simulation-specific information needed to run the full model. This document lists the available annotations for the various `SBML components <http://sbml.org/Documents/Specifications#SBML_Level_2>`_. All parameters which are given as `MathML <http://www.w3.org/Math/>`_ expressions can make use of the model :ref:`parameters <parameters>`.

SBML Container
--------------

* ``im:sbml``

  * ``im:solver`` The solver type. Available types are ``stochastic_homogeneous``, ``stochastic_inhomogeneous``, ``deterministic_homogeneous`` and ``deterministic_inhomogeneous``.
  * ``im:runs`` The number of runs as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:time`` The simulation run time as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:outputInterval`` The output interval as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:init`` The :ref:`initialization method <initialization>`.

    * ``im:script``  The actual script as a `XML CDATA array <http://www.w3schools.com/xml/xml_cdata.asp>`_.

  * ``im:events`` The :ref:`events method <events>`.

    * ``im:script`` The actual script as a `XML CDATA array <http://www.w3schools.com/xml/xml_cdata.asp>`_.
    * ``im:listOfTimeStamps`` A list of the time stamps for which the events script is called.

  * ``im:driftDiffusivity`` The :ref:`inhomogeneous drift-diffusivity method <drift_diffusivity_method>`. Attributes are ``nonLinear`` and ``computeMoments`` which can be set to ``true`` or ``false``.
 
    * ``im:method`` The actual method as a `XML CDATA array <http://www.w3schools.com/xml/xml_cdata.asp>`_.

  * ``im:setField`` The :ref:`update fields <update_fields>` method.
 
    * ``im:method`` The actual method as a `XML CDATA array <http://www.w3schools.com/xml/xml_cdata.asp>`_.

  * ``im:newIndividualsMethod`` The :ref:`New Individuals Method <new-individuals>`.
 
    * ``im:method`` The actual method as a `XML CDATA array <http://www.w3schools.com/xml/xml_cdata.asp>`_.

Model
-----

* ``im:model`` The Inchman model annotations.

  * ``im:gridWidth`` The grid width as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:gridHeight`` The grid height as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:physicalWidth`` The physical width as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:physicalHeight`` The physical height as a `MathML <http://www.w3.org/Math/>`_ expression.

Compartment
-----------

* ``im:compartment`` The Inchman compartment annotations.

  * ``im:x`` The compartment position (x) as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:y`` The compartment position (y) as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:width`` The compartment width as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:height`` The compartment height as a `MathML <http://www.w3.org/Math/>`_ expression.

Species
-------

.. note :: 

  Inchman does not support constraining species to particular compartments. All species have the obligatory ``World`` compartment as the associated compartment.

* ``im:species`` The Inchman species annotations.

  * ``im:diffusionConstant`` The diffusivity of the species as a `MathML <http://www.w3.org/Math/>`_ expression.
  * ``im:individual`` If this tag is present, the species is an individual species.
  * ``im:listOfCompartmentParameters`` Can be used to specify the initial amount of this species in the compartment

    * ``im:initialAmount`` The initial amount of this species in this compartment.

Parameters
----------

.. note :: 

  Parameters are defined using the SBML ``parameter`` tag with the corresponding ``id``. The ``value`` tag is ignored and instead the parameter value is set from the information specified in the annotation. The value of :ref:`field parameters <field_parameter>` is set in the :ref:`initialization method <initialization>`.

* ``im:parameter`` The Inchman parameter annotations.

  * ``type`` This attribute can be set to either ``int`` or ``float``
  * ``domain`` This attribute can be set to ``single``, ``range``, ``random`` or ``field``.
  * ``from`` The bottom end of the sweeping range. If the parameter domain is ``single``, this is the value of the parameter.
  * ``to`` The top end of the sweeping range.
  * ``step`` The step size for the sweeping.
  * ``points`` The number of sample points for the sweeping.

Reactions
---------

* ``im:reaction`` The Inchman reaction annotations.

  * ``im:listOfCompartmentParameters`` Can be used to constrain the reaction to particular compartments.

    * ``im:compartmentParameters`` The compartments this reaction is constrained to are given in the attribute ``compartment``
  
      * ``im:isAllowed`` A `MathML <http://www.w3.org/Math/>`_ expression evaluating to either ``true`` or ``false`` to allow or disallow this reaction in the compartment.
