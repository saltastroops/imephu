# Development Deployment

This job deploys the package to a development PyPI server.

There are two differences to the production deployment:

1. A username and password are used for authentication.
2. The server is a development server.

## Requirements

Before using the Jenkinsfile, you must create credentials in Jenkins for the PyPI server user. These must be of the type "Username and password", and their id must be "pypi-dev-credentials".

