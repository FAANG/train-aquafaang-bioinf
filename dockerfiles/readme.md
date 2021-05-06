# Pull Repo Container Images 
To pull the container images associated with the repository from Github Container Registry there are some requirements

- Github account.
- Authentication with a personal access token with read:packages scope. [Instructions for generating this token](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token#:~:text=Creating%20a%20token,-Verify%20your%20email&text=In%20the%20upper%2Dright%20corner,Click%20Generate%20new%20token.).

To pull the image: 
- Save your PAT as an environment variable.

    ```sh
    export CR_PAT=YOUR_TOKEN
    ```

- Sign in to the Container registry service at ghcr.io
    ```sh
    echo $CR_PAT | docker login ghcr.io -u GITHUB_USERNAME --password-stdin
     ```
    ```sh
    docker pull ghcr.io/OWNER/IMAGE_NAME
    ```
- For the container images associated with this repo:
    ```sh
    docker pull ghcr.io/faang/train-aquafaang-genrich:latest
    docker pull ghcr.io/faang/train-aquafaang-samtools:latest
    ```