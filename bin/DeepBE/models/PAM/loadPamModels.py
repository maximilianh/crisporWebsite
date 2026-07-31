import os
import tensorflow as tf

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def loadModel():

    model0 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_SpCas9_model.h5"), compile=False
    )
    model1 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_VRQR_model.h5"), compile=False
    )
    model2 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_NG_model.h5"), compile=False
    )
    model3 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_NRRH_model.h5"), compile=False
    )
    model4 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_NRTH_model.h5"), compile=False
    )
    model5 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_NRCH_model.h5"), compile=False
    )
    model6 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_SpG_model.h5"), compile=False
    )
    model7 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_SpRY_model.h5"), compile=False
    )
    model8 = tf.keras.models.load_model(
        os.path.join(SCRIPT_DIR, "PAM_variant_Sc++_model.h5"), compile=False
    )
    model_list = [
        model0,
        model1,
        model2,
        model3,
        model4,
        model5,
        model6,
        model7,
        model8,
        ]

    return model_list
