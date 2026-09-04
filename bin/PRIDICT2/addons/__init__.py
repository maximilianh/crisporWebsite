"""PRIDICT2.0 add-ons that generate batch input sequences.

The two sub-packages turn a single target sequence into an iterable of
PRIDICT2.0 batch inputs which can be passed directly to
``pridict2_pegRNA_design.predict_batch_sequences()``::

    from addons.flexible_mutations import flexible_mutation_sequences
    from addons.silentbystander import silent_bystander_sequences
"""
