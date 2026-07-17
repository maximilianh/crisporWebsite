- For new installations, two virtual environments needs to be installed to run Base editing scoring models.
- please install them by running the commands below


## Notes

- DeepBE and FORECast-BE were designed with python 3.6, but some wheels aren't available anymore in this version.
- Here, the program is installed with python 3.9.25.
- After some testing, the outputs appear to be identical using both versions of python.

# virtual evironment for DeepBE

- see https://github.com/NahyeKim/DeepBE
   
``` 
    python3.9 -m venv venvDeepBE
    pip install pandas==1.3.0 tensorflow==2.6.2 protobuf==3.20.3 tensorRT
```

# virtual environment for FORECasT-BE

- see https://github.com/ananth-pallaseni/FORECasT-BE
- FORECasT-BE depends on scikit-learn 0.23, but it isn't availble for python3.9, so version 0.24.1 is used instead
- Ather some testing, both versions produce the same results

```
    python3.9 -m venv venvForecastBe
    pip install biopython==1.85 numpy==1.23.5 pandas==2.3.3 scikit-learn==0.24.1 scipy==1.13.1
```
