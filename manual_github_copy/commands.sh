pip uninstall leap_gwas

python setup.py sdist upload -r pypitest
pip install --user  -i https://testpypi.python.org/pypi leap_gwas

#python setup.py sdist upload -r pypi
#pip install --user  -i https://pypi.python.org/pypi leap_gwas