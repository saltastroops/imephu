# Setting up the Jenkins CI/CD workflow

## Tunneling through a firewall

In case your Jenkins server is behind a firewall, you may use [ngrok](https://ngrok.io) to allow the use of GitHub webhooks. See [this video](https://youtu.be/yMNJeWeE0qI) for a walk-through. ngrok lets you create a permanent domain, and you can specify it when creating the tunnel:

```bash
ngrok http --domain your.permanent.domain 8080
```

`your.permanent.domain` and `8080` need to be replaced with your domain and your Jenkins server's port, respectively.

## Prerequisites

### GitHub App

Before setting up the Jenkins job you need to create a GitHub App, install it in your repository, and add credentials (of type "GitHub App") for it in Jenkins. An explanation is given in [this video](https://youtu.be/aDmeeVDrp0o) and on [this webpage](https://docs.cloudbees.com/docs/cloudbees-ci/latest/cloud-admin-guide/github-app-auth). When creating the app, you are asked to generate a private key. This key will automatically be downloaded. Note that, as mentioned in these resources, you need to convert the downloaded key with a command like

```bash
openssl pkcs8 -topk8 -inform PEM -outform PEM -in downloaded-github-app-key.pem -out converted-github-app-key.pem -nocrypt
```

You need to create credentials of type "GitHub App" for this GitHub App. The App ID is an integer, which you can find on the General tab of the app's settings. The Key is the converted private key you generated above.

![GitHub App credentials](img/github_app_credentials.png)

### Environment variables

The following environment variable needs to be defined.

| Environment variable | Description                                                                                   | Example value        |
|----------------------|-----------------------------------------------------------------------------------------------|----------------------|
| SALT_ASTROOPS_EMAIL  | Email address for SALT Astronomy Operations. This is used for pipeline failure notifications. | astroops@example.com |

Environment variables can be set via Manage Jenkins - Configure System (cf. [this Stack Overflow entry](https://stackoverflow.com/questions/54207815/does-jenkins-have-a-feature-like-credentials-for-non-secrets)).

### Credentials

The following credentials need to be defined.

| Credentials id          | Credentials type       | Description                                                                                   |
|-------------------------|------------------------|-----------------------------------------------------------------------------------------------|
| saltastroops_dev_pypi   | Username with password | Credentials for authenticating on the development PyPI server (https://pypi.cape.saao.ac.za). | 
| saltastroops_pypi_token | Secret string          | Token for authenticating the user `saltastroops` on the PyPI server (https://pypi.org).       | 

### Email Extension plugin

The Email Extension plugin must be installed, and it must have been configured in Manage Jenkins - System. See [this webpage](https://www.edureka.co/blog/email-notification-in-jenkins/) for more details.

## Setting up the workflow

Go to the Jenkins dashboard and click on "New item" in the sidebar menu. Choose a name and select "Multibranch Pipeline" as the type of item to create.

On the configuration page add a GitHub source.

![Add a GitHub source](img/add-github-source.png)

Choose the GitHub App credentials defined above as the credentials and https://github.com/saltastroops/imephu as the repository HTTPS URL.

![Set the credentials and repository](img/github-credentials-and-repo.png)

Still for the GitHub source, add the following behaviours by using the "Add" button:

* Filter by name (with wildcards)
* Clean after checkout
* Clean before checkout

For the filter enter `main development` in the Include input and leave the Exclude input empty.

![Add the behaviours](img/github-source-behaviours.png)

In the Build Configuration section enter `jenkins/Jenkinsfile` as the script path.

![Choose the build configuration path](img/build-configuration.png)

Finally, if webhooks are unavailable (either because they are not enabled on GitHub or because the Jenkins server cannot be connected to from the outside), you should enable periodic scanning of the repository.

![Enable periodic scanning of the repository](img/polling.png)

## Updating the workflow when supported Python versions change

When the supported Python versions change, you have to update two files.

1. In `jenkins/Dockerfile` you have to update the Python versions that are installed with pyenv.
2. In `jenkins/Jenkinsfile` you have to update the Python versions in the `versions` array and (if necessary) the value of the environment variable `PYTHON_VERSION`, which is defined in the `environment` directive.

The Python versions in the Dockerfile and Jenkinsfile must be identical.

The Python versions installed with pyenv are relevant for running Python; the linting checks are run with the system Python version defined by the Docker base image.
