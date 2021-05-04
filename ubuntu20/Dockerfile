FROM ubuntu:20.10

# Env & Arg variables
ARG USERNAME=user
ARG USERPASS=mypassword

# Apt update & apt install required packages
# whois: required for mkpasswd
RUN apt-get update && apt-get -y install openssh-server whois

# Add a non-root user & set password
RUN useradd -ms /bin/bash $USERNAME

# Set password for non-root user
RUN usermod --password $(echo "$USERPASS" | mkpasswd -s) $USERNAME

# Create the ssh directory and authorized_keys file
USER $USERNAME
RUN mkdir /home/$USERNAME/.ssh && touch /home/$USERNAME/.ssh/authorized_keys
USER root

# ------------------------------------------------------------------------------------------------------
RUN apt-get install -y build-essential gfortran libnetcdf-dev libgsl-dev python2 gawk python3-pip
# RUN apt-get install -y libgeos-dev libproj-dev liblapack-dev libblas-dev gfortran
# update again to ensure
RUN apt-get update && apt-get upgrade
# ------------------------------------------------------------------------------------------------------

# Remove no-needed packages
RUN apt-get purge -y whois && apt-get -y autoremove && apt-get -y autoclean && apt-get -y clean

# Set volumes
VOLUME /home/$USERNAME/.ssh
VOLUME /etc/ssh

# expose ports
EXPOSE 22

# Copy the entrypoint
COPY ./entrypoint.sh /root/entrypoint.sh
RUN chmod +x /root/entrypoint.sh

# Run entrypoint
CMD ["/root/entrypoint.sh"]