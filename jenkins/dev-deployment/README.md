# Development Deployment

This job deploys the package to a development PyPI server.

There are two differences to the production deployment:

1. A username and password are used for authentication.
2. The server is not hardcoded but must be passed as a parameter.

## Requirements

Before using the Jenkinsfile, you must do the following.

1. Create credentials in Jenkins for the PyPI server user. These must be of the type "Username and password".
2. When 
