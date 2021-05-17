from setuptools import setup

setup(
  name='scGene2Vec',
  version='0.1.0',
  author='Beth Baumann',
  author_email='baumann.bethany@gmail.com',
  packages=['scGene2Vec'],
  description='trying single cell gene vectors with gensim',
  long_description=open('README.txt').read(),
  install_requires=[
     "pandas >= 1.2.2",
     "gensim >= 3.8.0",
 ],
)
