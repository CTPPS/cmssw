from weakref import WeakKeyDictionary
import automation_control as ctrl
import argparse
import enum 
import logging
from typing import Any, Type, Union
from os import listdir, walk, environ
from os.path import isfile, join

logger = logging.getLogger("EfficiencyAnalysisLogger")
logger.setLevel(logging.DEBUG)

ch = logging.FileHandler("EfficiencyAnalysisEngine.log")
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

campaign=environ.get("CAMPAIGN")
workflow=environ.get("WORKFLOW")
dataset=environ.get("DATASET")
proxy=environ.get("PROXY")

template_for_first_module = "CrabConfigTemplateForFirstModule.py"

@ctrl.define_status_enum
class TaskStatusEnum(enum.Enum):
    """
        Class to encode enum tasks statuses for the purpouse of this automation workflow
    """
    initialized = enum.auto(),
    duringFirstWorker = enum.auto(),
    waitingForFirstWorkerTransfer= enum.auto()
    done = enum.auto()
    

@ctrl.decorate_with_enum(TaskStatusEnum)
class TaskStatus:
    loop_id = 0.0    
    condor_job_id = 0

def get_tasks_numbers_list(tasks_list_path):
    with open(tasks_list_path) as tasks_list_path:
        tasks_list_data = tasks_list_path.read()
        tasks_list_data = tasks_list_data.replace(" ", "")
        tasks_list = tasks_list_data.split(",")
    return tasks_list


def prepare_parser()->argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=
    """This is a script to run PPS Efficiency Analysis automation workflow""", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-t', '--tasks_list', dest='tasks_list_path', help='path to file containing list of data periods', required=True)   
    return parser


def get_runs_range(data_period):
    """MOCKED"""
    return '317080'


def process_new_tasks(tasks_list_path, task_controller):
    tasks_list = get_tasks_numbers_list(tasks_list_path)
    tasks_list = set(tasks_list)
    tasks_in_database = task_controller.getAllTasks().get_points()
    tasks_in_database = set(map(lambda x: x['dataPeriod'], tasks_in_database))
    tasks_not_submited_yet = tasks_list-tasks_in_database
    if tasks_not_submited_yet:
        task_controller.submitTasks(tasks_not_submited_yet)
   

def submit_task_to_crab(campaign, workflow, data_period, dataset, template, proxy):
    result = ctrl.submit_task_to_crab(campaign, workflow, data_period, get_runs_range(data_period), template, dataset, proxy)

    return result


def set_status_after_first_worker_submission(task_status, operation_result):
    task_status.duringFirstWorker=1
    task_status.initialized=0
    task_status.loop_id+=1
    return task_status

             
storage_path = "/eos/user/m/mobrzut"

TRANSITIONS_DICT = {
                        'initialized': (submit_task_to_crab, 0, set_status_after_first_worker_submission, [dataset, template_for_first_module, proxy] ),
                        'duringFirstWorker': (ctrl.check_if_crab_task_is_finished, True, TaskStatus.waitingForFirstWorkerTransfer, [proxy]),
                        'waitingForFirstWorkerTransfer': (ctrl.is_crab_output_already_transfered, True, TaskStatus.done, [proxy])
                   } 


    
if __name__ == '__main__':
    parser = prepare_parser()
    opts = parser.parse_args()
    task_controller = ctrl.TaskCtrl.TaskControl(campaign=campaign, workflow=workflow, TaskStatusClass=TaskStatus)
    process_new_tasks(opts.tasks_list_path, task_controller)
    finite_state_machine = ctrl.FiniteStateMachine(TRANSITIONS_DICT)
    finite_state_machine.process_tasks(task_controller, TaskStatusClass=TaskStatus)
    
    
    