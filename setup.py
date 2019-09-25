import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name='molreps',
        version="0.0.1",
        description='Molecular representation for machine learning',
        author='Dakota Folmsbee',
        author_email='dfolmsbee@gmail.com',
        url="https://github.com/dlf57/MolecularLearning",
        packages=setuptools.find_packages(),
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3',
        ],
    )