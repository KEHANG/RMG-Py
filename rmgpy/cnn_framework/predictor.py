
from cnn_model import build_model
from .input import read_input_file
import os
import rmgpy

class Predictor(object):

	def __init__(self, input_file=None):

		self.model = None
		if input_file:
			self.input_file = input_file
		else:
			self.input_file = os.path.join(os.path.dirname(rmgpy.__file__),
										'cnn_framework',
										'test_data', 
										'minimal_predictor', 
										'predictor_input.py'
										)

	def build_model(self):

		self.model = build_model()

	def load_input(self, path=None):
		
		if path is None: 
			path = self.input_file
			print path
		read_input_file(path, self)

	def train(self):

		pass

	def load_parameters(self):

		pass

	def save_parameters(self):

		pass

	def predict(self):

		pass
    
