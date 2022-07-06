from django.db import models
from django.urls import reverse
from .mol_validate import MolValidate

# Create your models here.

class Molecule(models.Model):
	"""Model representing a single molecule that someone has done a search on"""
	SMI = models.CharField(max_length=500)
	smiID = models.CharField(max_length=100)

	def __str__(self):
		return self.smiID
