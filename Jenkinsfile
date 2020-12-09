pipeline {
  agent any
  stages {
    stage('Checkout') {
      steps {
        catchError(buildResult: 'FAILURE', message: 'Failed at the checkout step', stageResult: 'FAILURE') {
          git 'https://github.com/larson-group/clubb.git'
        }

      }
    }

    stage('Compile') {
      steps {
        dir(path: 'clubb/compile/') {
          catchError(buildResult: 'FAILURE', message: 'Failed on the compile step', stageResult: 'FAILURE') {
            sh './compile.bash -c config/linux_x86_64_gfortran.bash'
          }

        }

      }
    }

    stage('Run') {
      steps {
        catchError(buildResult: 'FAILURE', message: 'Failed to run clubb', stageResult: 'FAILURE') {
          dir(path: 'clubb/run_scripts/') {
            sh './run_scm_all.bash'
          }

        }

      }
    }

  }
}