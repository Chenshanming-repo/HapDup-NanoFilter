import sys
import torch
import os
import h5py
from os.path import isfile, join
from os import listdir
import time
import torch.distributed as dist
import torch.nn as nn
import torch.multiprocessing as mp
from datetime import datetime
# Custom generator for our dataset
from torch.utils.data import DataLoader
from pepper_variant.modules.python.models.dataloader import SequenceDataset
from pepper_variant.modules.python.models.ModelHander import ModelHandler
from pepper_variant.modules.python.models.test import test
from pepper_variant.modules.python.Options import ImageSizeOptions, ImageSizeOptionsHP, TrainOptions

os.environ['PYTHONWARNINGS'] = 'ignore:semaphore_tracker:UserWarning'

"""
Train a model and return the model and optimizer trained.

Input:
- A train CSV containing training image set information (usually chr1-18)

Return:
- A trained model
"""


def save_best_model(transducer_model, model_optimizer, hidden_size, layers, epoch,
                    file_name):
    """
    Save the best model
    :param transducer_model: A trained model
    :param model_optimizer: Model optimizer
    :param hidden_size: Number of hidden layers
    :param layers: Number of GRU layers to use
    :param epoch: Epoch/iteration number
    :param file_name: Output file name
    :return:
    """
    if os.path.isfile(file_name):
        os.remove(file_name)
    ModelHandler.save_checkpoint({
        'model_state_dict': transducer_model.state_dict(),
        'model_optimizer': model_optimizer.state_dict(),
        'hidden_size': hidden_size,
        'gru_layers': layers,
        'epochs': epoch,
    }, file_name)
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL" + file_name + " SAVED SUCCESSFULLY.\n")
    sys.stderr.flush()


def train(train_file, test_file, batch_size, test_batch_size, step_size, epoch_limit, gpu_mode, num_workers, retrain_model,
          retrain_model_path, gru_layers, hidden_size, lr, decay, model_dir, stats_dir, hp_mode, train_mode):
    train_loss_logger = open(stats_dir + "train_loss.csv", 'w')
    test_loss_logger = open(stats_dir + "test_loss.csv", 'w')
    confusion_matrix_logger = open(stats_dir + "base_confusion_matrix.txt", 'w')

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: LOADING DATA\n")
    sys.stderr.flush()

    train_with_base = False

    train_data_set = SequenceDataset(train_file)

    # shuffle is off because we are using distributed sampler, which has shuffle on
    train_loader = torch.utils.data.DataLoader(
        dataset=train_data_set,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True)

    if not hp_mode:
        image_height = ImageSizeOptions.IMAGE_HEIGHT
        num_classes = ImageSizeOptions.TOTAL_LABELS
        num_type_classes = ImageSizeOptions.TOTAL_TYPE_LABELS
    else:
        image_height = ImageSizeOptionsHP.IMAGE_HEIGHT
        num_classes = ImageSizeOptionsHP.TOTAL_LABELS
        num_type_classes = ImageSizeOptionsHP.TOTAL_TYPE_LABELS

    if retrain_model is True:
        if os.path.isfile(retrain_model_path) is False:
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: INVALID PATH TO RETRAIN PATH MODEL --retrain_model_path\n")
            exit(1)
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RETRAIN MODEL LOADING\n")
        transducer_model, hidden_size, gru_layers, prev_ite = \
            ModelHandler.load_simple_model_for_training(retrain_model,
                                                        image_features=image_height,
                                                        num_classes=num_classes,
                                                        num_type_classes=num_type_classes)

        if train_mode is True:
            epoch_limit = prev_ite + epoch_limit

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RETRAIN MODEL LOADED\n")
    else:
        transducer_model = ModelHandler.get_new_gru_model(image_features=image_height,
                                                          gru_layers=gru_layers,
                                                          hidden_size=hidden_size,
                                                          num_classes=num_classes,
                                                          num_classes_type=num_type_classes)
        prev_ite = 0

    param_count = sum(p.numel() for p in transducer_model.parameters() if p.requires_grad)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL TRAINABLE PARAMETERS:\t" + str(param_count) + "\n")
    sys.stderr.flush()
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL OPTIMIZER PARAMETERS:\n")
    sys.stderr.flush()
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: LEARNING RATE: " + str(lr) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: WEIGHT DECAY: " + str(decay) + "\n")

    model_optimizer = torch.optim.Adam(transducer_model.parameters(), lr=lr, weight_decay=decay)
    # scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(model_optimizer, 'min')

    if retrain_model is True:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: OPTIMIZER LOADING\n")
        model_optimizer = ModelHandler.load_simple_optimizer(model_optimizer, retrain_model_path, gpu_mode)
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: OPTIMIZER LOADED\n")
        sys.stderr.flush()

    if gpu_mode:
        transducer_model = transducer_model.cuda()
        transducer_model = nn.parallel.DataParallel(transducer_model)

    weights = [0.3, 1.0, 1.0]
    class_weights = torch.tensor(weights).cuda()
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: Class weights: " + str(weights) + "\n")
    # Loss
    if train_with_base:
        criterion_base = nn.CrossEntropyLoss()
    criterion_type = nn.CrossEntropyLoss(weight=class_weights)

    if gpu_mode is True:
        if train_with_base:
            criterion_base = criterion_base.cuda()
        criterion_type = criterion_type.cuda()

    start_epoch = prev_ite

    # Train the Model
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TRAINING STARTING\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: START: " + str(start_epoch + 1) + " END: " + str(epoch_limit) + "\n")
    sys.stderr.flush()

    # Creates a GradScaler once at the beginning of training.
    # scaler = torch.cuda.amp.GradScaler()

    start_iteration = 0
    if step_size == 0:
        step_size = len(train_loader)

    iteration_limit = int((step_size * epoch_limit) / len(train_loader)) + 1
    step_no = 0
    epoch = 0

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ITERATION LIMIT: " + str(iteration_limit - 1) + "\n")
    sys.stderr.flush()

    for iteration in range(start_iteration, iteration_limit, 1):
        start_time = time.time()
        total_loss = 0
        total_base_loss = 0
        total_type_loss = 0
        total_images = 0

        for images, base_labels, type_labels in train_loader:
            # make sure the model is in train mode.
            transducer_model = transducer_model.train()

            # image and label handling
            if train_with_base:
                base_labels = base_labels.type(torch.LongTensor)
            type_labels = type_labels.type(torch.LongTensor)
            images = images.type(torch.FloatTensor)
            if gpu_mode:
                images = images.cuda()
                if train_with_base:
                    base_labels = base_labels.cuda()
                type_labels = type_labels.cuda()

            # generate hidden and cell state
            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
            cell_state = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            if gpu_mode:
                hidden = hidden.cuda()
                cell_state = cell_state.cuda()

            # set the optimizer to zero grad
            model_optimizer.zero_grad()

            output_base, output_type = transducer_model(images, hidden, cell_state, train_mode)

            if train_with_base:
                loss_base = criterion_base(output_base.contiguous().view(-1, num_classes), base_labels.contiguous().view(-1))
            loss_type = criterion_type(output_type.contiguous().view(-1, num_type_classes), type_labels.contiguous().view(-1))

            if train_with_base:
                base_type_loss = 0.3 * loss_base + loss_type
            else:
                base_type_loss = loss_type

            base_type_loss.backward()

            model_optimizer.step()

            # done calculating the loss
            total_loss += base_type_loss.item()
            if train_with_base:
                total_base_loss += loss_base.item()
            total_type_loss += loss_type.item()
            total_images += images.size(0)

            # update the progress bar
            avg_loss = (total_loss / total_images) if total_images else 0
            if train_with_base:
                avg_base_loss = (total_base_loss / total_images) if total_images else 0
            avg_type_loss = (total_type_loss / total_images) if total_images else 0

            if step_no % 100 == 0 and step_no > 0:
                nom = step_no - ((iteration * len(train_loader)) - (iteration * len(train_loader)) % 100)
                denom = len(train_loader)
                percent_complete = int((float(nom) / float(denom) * 100))
                time_now = time.time()
                mins = int((time_now - start_time) / 60)
                secs = int((time_now - start_time)) % 60

                if train_with_base:
                    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: "
                                     + " ITERATION: " + str(iteration + 1)
                                     + " EPOCH: " + str(epoch + 1)
                                     + " STEP: " + str(step_no) + "/" + str((epoch + 1) * step_size) + "/" + str((iteration + 1) * len(train_loader))
                                     + " LOSS: " + "{:.9f}".format(avg_loss)
                                     + " BASE LOSS: " + "{:.9f}".format(avg_base_loss)
                                     + " TYPE LOSS: " + "{:.9f}".format(avg_type_loss)
                                     + " EPOCH COMPLETE: (" + str(percent_complete) + "%)"
                                     + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                else:
                    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: "
                                     + " ITERATION: " + str(iteration + 1)
                                     + " EPOCH: " + str(epoch + 1)
                                     + " STEP: " + str(step_no) + "/" + str((epoch + 1) * step_size) + "/" + str((iteration + 1) * len(train_loader))
                                     + " LOSS: " + "{:.9f}".format(avg_loss)
                                     + " EPOCH COMPLETE: (" + str(percent_complete) + "%)"
                                     + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                sys.stderr.flush()

            if step_no % step_size == 0 and step_no > 0:
                transducer_model = transducer_model.eval()
                torch.cuda.empty_cache()

                stats_dictionary = test(test_file, test_batch_size, gpu_mode, transducer_model, num_workers,
                                        gru_layers, hidden_size, train_with_base, num_classes=ImageSizeOptions.TOTAL_LABELS)

                save_best_model(transducer_model, model_optimizer, hidden_size, gru_layers, epoch, model_dir + "PEPPER_VARIANT_STEP_" + str(step_no) + '_checkpoint.pkl')
                train_loss_logger.write(str(step_no) + "," + str(avg_loss) + "\n")

                test_loss_logger.write(str(step_no) + "," + str(stats_dictionary['loss']) + "," + str(stats_dictionary['base_accuracy']) + "\n")
                confusion_matrix_logger.write(str(step_no) + "\n")

                # print confusion matrix to file
                if train_with_base:
                    confusion_matrix_logger.write("Base confusion Matrix:" + "\n")
                    confusion_matrix_logger.write("            ")
                    for label in ImageSizeOptions.decoded_base_labels:
                        confusion_matrix_logger.write(str(label) + '         ')
                    confusion_matrix_logger.write("\n")

                    for i, row in enumerate(stats_dictionary['base_confusion_matrix'].value()):
                        confusion_matrix_logger.write(str(ImageSizeOptions.decoded_base_labels[i]) + '   ')
                        for j, val in enumerate(row):
                            confusion_matrix_logger.write("{0:9d}".format(val) + '  ')
                        confusion_matrix_logger.write("\n")
                    confusion_matrix_logger.flush()

                confusion_matrix_logger.write("Type confusion Matrix:" + "\n")
                confusion_matrix_logger.write("            ")
                for label in ImageSizeOptions.decoded_labels:
                    confusion_matrix_logger.write(str(label) + '         ')
                confusion_matrix_logger.write("\n")

                for i, row in enumerate(stats_dictionary['type_confusion_matrix'].value()):
                    confusion_matrix_logger.write(str(ImageSizeOptions.decoded_base_labels[i]) + '   ')
                    for j, val in enumerate(row):
                        confusion_matrix_logger.write("{0:9d}".format(val) + '  ')
                    confusion_matrix_logger.write("\n")
                confusion_matrix_logger.flush()

                train_loss_logger.flush()
                test_loss_logger.flush()
                confusion_matrix_logger.flush()

                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST COMPLETED.\n")
                sys.stderr.flush()
                transducer_model = transducer_model.train()

                epoch += 1
                # scheduler.step(stats_dictionary['loss'])

            # increase step size
            step_no += 1

            if epoch == epoch_limit:
                break

        time_now = time.time()
        mins = int((time_now - start_time) / 60)
        secs = int((time_now - start_time)) % 60
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ELAPSED TIME FOR ONE ITERATION: " + str(mins) + " Min " + str(secs) + " Sec\n")

        if epoch == epoch_limit:
            break

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED TRAINING\n")


def cleanup():
    dist.destroy_process_group()


def setup(rank, device_ids, args):
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '12355'

    # initialize the process group
    dist.init_process_group("gloo", rank=rank, world_size=len(device_ids))

    train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model, \
    retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir, stats_dir, total_callers, \
    train_mode = args

    # issue with semaphore lock: https://github.com/pytorch/pytorch/issues/2517
    # mp.set_start_method('spawn')

    # Explicitly setting seed to make sure that models created in two processes
    # start from same random weights and biases. https://github.com/pytorch/pytorch/issues/2517
    torch.manual_seed(42)
    train(train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model, retrain_model_path,
          gru_layers, hidden_size, learning_rate, weight_decay, model_dir, stats_dir, train_mode,
          total_callers, rank, device_ids[rank])
    cleanup()

def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-3:] == 'hdf']
    return file_paths


def train_distributed(train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model,
                      retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir,
                      stats_dir, device_ids, total_callers, hp_mode, train_mode):

    # run linearly
    train(train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model, retrain_model_path,
          gru_layers, hidden_size, learning_rate, weight_decay, model_dir, stats_dir, hp_mode, train_mode)

    # args = (train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model,
    #         retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir,
    #         stats_dir, total_callers, train_mode)
    #
    # mp.spawn(setup,
    #          args=(device_ids, args),
    #          nprocs=len(device_ids),
    #          join=True)
