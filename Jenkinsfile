/*
###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
###############################################################################

This file defines the build/test/post steps for Jenkins.

Dependencies: failed-test-email.template  // This parses the output of Jenkins and formats the email that Jenkins sends
*/

def email_to = "idaes.jenkins@lbl.gov"  // The email address that the build email will go to
def email_reply_to = "mrshepherd@lbl.gov"  // The email address that will be in the reply-to

pipeline {
  agent { 
    docker { 
      image 'conda/miniconda3-centos7:latest'
    } 
  }
  stages {
    // Commented out for an example until we start using these parameters
    // stage('cron-nightly-test') {
    //   when {
    //     expression { params.BUILD_SCHEDULE == 'Nightly'}
    //   }
    //   steps {
    //     sh 'echo "nightly works"'
    //   }
    // }
    // stage('cron-weekly-test') {
    //   when {
    //     expression { params.BUILD_SCHEDULE == 'Weekly'}
    //   }
    //   steps {
    //     sh 'echo "weekly works"'
    //   }
    // }
    stage('root-setup') {
      steps {
        // slackSend (message: "Build Started - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
        sh 'mkdir -p $JENKINS_HOME/email-templates'
        sh 'cp failed-test-email.template $JENKINS_HOME/email-templates'
        sh 'yum install -y gcc g++ git gcc-gfortran libboost-dev make'
        sh 'pwd'

      }
    }
    stage('3.6-setup') {
      when {
        expression { params.PYTHON_VERSION == '3.6'}
      }
      steps {
        sh '''
         conda create -n idaes3.6 python=3.6 pytest
         source activate idaes3.6
         pip install -r requirements-dev.txt --user jenkins
         export TEMP_LC_ALL=$LC_ALL
         export TEMP_LANG=$LANG
         export LC_ALL=en_US.utf-8
         export LANG=en_US.utf-8
         python setup.py install
         idaes get-extensions
         export LC_ALL=$TEMP_LC_ALL
         export LANG=$TEMP_LANG
         source deactivate
         '''
      }
    }
    stage('3.7-setup') {
      when {
        expression { params.PYTHON_VERSION == '3.7'}
      }
      steps {
        sh '''
         conda create -n idaes3.7 python=3.7 pytest
         source activate idaes3.7
         pip install -r requirements-dev.txt --user jenkins
         export TEMP_LC_ALL=$LC_ALL
         export TEMP_LANG=$LANG
         export LC_ALL=en_US.utf-8
         export LANG=en_US.utf-8
         python setup.py install
         idaes get-extensions
         export LC_ALL=$TEMP_LC_ALL
         export LANG=$TEMP_LANG
         source deactivate
         '''
      }
    }
    stage('3.6-test') {
      when {
        expression { params.PYTHON_VERSION == '3.6'}
      }
      steps {
        catchError(buildResult: 'FAILURE', stageResult: 'FAILURE') {
          sh '''
           source activate idaes3.6
           pylint -E --ignore-patterns="test_.*" idaes || true
           pytest --junitxml=results.xml -c pytest.ini idaes
           source deactivate
           '''
        }   
      }
    }
    stage('3.7-test') {
      when {
        expression { params.PYTHON_VERSION == '3.7'}
      }
      steps {
        catchError(buildResult: 'FAILURE', stageResult: 'FAILURE') {
          sh '''
           source activate idaes3.7
           pylint -E --ignore-patterns="test_.*" idaes || true
           pytest --junitxml=results.xml -c pytest.ini idaes
           source deactivate
           '''
        }
      }   
    }
  }
  post {
    always {
      junit allowEmptyResults: true, testResults: 'results.xml'
      emailext attachLog: true, body: '''${SCRIPT, template="failed-test-email.template"}''', replyTo: '${email_reply_to}',
       subject: "${JOB_NAME} ${params.PYTHON_VERSION} - Build ${BUILD_NUMBER} ${currentBuild.result}", to: '${email_to}'
    }
    // success {
    //   slackSend (color: '#00FF00', message: "SUCCESSFUL - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
    //   emailext attachLog: true, body: "${currentBuild.result}: ${BUILD_URL}", compressLog: true, replyTo: 'mrshepherd@lbl.gov',
    //    subject: "Build Log: ${JOB_NAME} - Build ${BUILD_NUMBER} ${currentBuild.result}", to: 'mrshepherd@lbl.gov'
    // }

    // failure {
    //   slackSend (color: '#FF0000', message: "FAILED - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
    //   emailext attachLog: true, body: "${currentBuild.result}: ${env.BUILD_URL}", compressLog: true, replyTo: 'mrshepherd@lbl.gov',
    //    subject: "Build Log: ${env.JOB_NAME} - Build ${env.BUILD_NUMBER} ${currentBuild.result}", to: 'mrshepherd@lbl.gov'
    // }
  }
}
