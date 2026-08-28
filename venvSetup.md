- For new installations, two virtual environments needs to be installed to run Base editing scoring models.
- please install them by running the commands below


## Notes

- DeepBE and FORECast-BE were designed with python 3.6, but some wheels aren't available anymore in this version.
- Here, the program is installed with python 3.9.25.
- After testing, the outputs appear to be identical using both versions of python.

# virtual evironment for DeepBE

- see https://github.com/NahyeKim/DeepBE
   
``` 
    python3.9 -m venv venvDeepBE
    source venvDeepBE/bin/activate
    pip install pandas==1.3.0 tensorflow==2.6.2 protobuf==3.20.3 tensorRT
    deactivate
```

# virtual environment for FORECasT-BE

- see https://github.com/ananth-pallaseni/FORECasT-BE
- FORECasT-BE depends on scikit-learn 0.23, but it isn't availble for python3.9, so version 0.24.1 is used instead
- Ather testing, both versions produce the same results

```
    python3.9 -m venv venvForecastBe
    source venvForecastBe/bin/activate
    pip install biopython==1.85 numpy==1.23.5 pandas==2.3.3 scikit-learn==0.24.1 scipy==1.13.1
    deactivate
```

# virtual environment for PRIDICT2.0

- see https://github.com/uzh-dqbm-cmi/PRIDICT2
- PRIDICT2.0 depends on python 3.10.12, but runs correctly on python 3.9.25 (TO VERRIFY)

 ```
    python3.9 -m venv venvPRIDICT2
    source venvPRIDICT2/bin/activate
    pip install tensorflow==2.13.1 biopython==1.81 joblib==1.3.1 matplotlib==3.7.2 pandas==2.0.3 prettytable==3.8.0 primer3-py==2.3.0 scikit-learn==1.3.0 scipy==1.11.1 seaborn==0.12.2 tqdm==4.65.0
    pip install torch==2.0.1 --index-url https://download.pytorch.org/whl/cpu
    deactivate
```

