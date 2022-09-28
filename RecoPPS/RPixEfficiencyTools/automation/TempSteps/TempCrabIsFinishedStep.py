####################################################

# Boiler plate code for importing automation_control form parent directory :SCREAM: 
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

######################################################

from urllib3 import proxy_from_url
import automation_control as ctrl

number = 37

campaign = f'ctrl_camp{number}'
workflow = f'crtl_workflow{number}'
data_period = 'foo_36' # todo use correct dataperiod number 
proxy = '/afs/cern.ch/user/m/mobrzut/public/jenkins_proxy/jenkins.pem'


ctrl.check_if_crab_task_is_finished(campaign, workflow, data_period, proxy)