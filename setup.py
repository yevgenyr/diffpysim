import setuptools

def main():
    ######################
    pkg_name = 'diffpysim'
    alias = "diffpysim"
    rel_ex_path = [pkg_name, "main"]
    version = "1.0_public"
    url = "git@github.com:yevgenyr/diffpysim.git"
    ######################

    # read long description from README
    with open("README.md", "r") as fh:
        long_description = fh.read()
    # create setuptools dict
    _setup = dict(
        name=pkg_name,
        version=version,
        author="Yevgeny Rakita",
        author_email="rakita@bgu.ac.il",
        description="diffpysim - simulating series of PDF and XRD",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url=url,
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: TBD",
            "Operating System :: OS Independent",
        ],
        entry_points={
            # define console_scripts here, see setuptools docs for details.
            'console_scripts': [
                f'{alias} = {".".join(rel_ex_path)}:main',
            ]
        },
        data_files=[("", ["LICENSE.txt", "README.md"])],
        include_package_data=False,
        # package_data={"": []},
        python_requires='>=3.7',
        zip_safe=False, )

    setuptools.setup(**_setup)


if __name__ == '__main__':
    main()
