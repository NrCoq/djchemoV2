from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import ugettext_lazy as _
from .mol_validate import MolValidate


class Search(forms.Form):
	smi = forms.CharField(label="", max_length=501, widget=forms.Textarea)

	def clean_smi(self):
		data = self.cleaned_data['smi']

		# Check if a date is not in the past.
		m = MolValidate
		if not m.valid_smi(data):
			raise ValidationError(_('Строка не имеет вид нотации SMILES'))

		mol = m.get_mol(data)
		smiles = m.get_std_smiles(mol)

		return smiles

