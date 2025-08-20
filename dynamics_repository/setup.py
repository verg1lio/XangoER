from setuptools import setup, find_packages

setup(
    name = "DynamicsSoftware",
    version = "0.1",
    packages=find_packages(),
    include_package_data=True,
    package_data ={"":["Software/Interface/saved_models.xlsx"]},
    install_requires=[
        "numpy",
        'matplotlib',
        'openpyxl',
        "pandas",
        "flet[all]"
    ]
)
