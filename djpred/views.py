from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from . import models
import os
from .mol_validate import MolValidate
from .forms import Search
from . import predictor
from . import similarity

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))



def landing(request):
	return render(request, 'landing.html')

@login_required
def index(request):

	if request.method == 'POST':
		form = Search(request.POST)

		if form.is_valid():
			data = form.cleaned_data['smi']
			m = MolValidate
			mol = m.get_mol(data)
			request.session['session_smiles'] = m.get_std_smiles(mol)
			return redirect('result')
		else:
			form = Search(request.POST)
			context = {
				'form':form
			}
			return render(request, 'index.html', context)
	elif request.method == 'GET':
		form = Search(request.GET)

		context = {
			'form': form
		}

		return render(request, 'index.html', context=context)

@login_required
def result(request):
	if 'session_smiles' not in request.session:
		return redirect('index')
	else:
		smiles = request.session['session_smiles']
		results = predictor.predict(smiles)

		sim = similarity.calc_similarity(smiles)

		intermed = Chem.MolFromSmiles(smiles)
		AllChem.Compute2DCoords(intermed)

		id = MolValidate.get_unique_id(intermed)

		Draw.MolToFile(intermed, os.path.join(BASE_DIR, 'djpred/static/' + id + '.svg'), size=(500, 500))
		pic_url = os.path.join('/static/' + id + '.svg')

		context = {
			'smiles': smiles,
			'nr_ahr': results['NR-AhR'],
			'nr_ar': results['NR-AR'],
			'nr_ar_lbd': results["NR-AR-LBD"],
			'nr_aromatase': results["NR-Aromatase"],
			'nr_er': results["NR-ER"],
			'nr_er_lbd': results["NR-ER-LBD"],
			'nr_ppar_gamma': results["NR-PPAR-gamma"],
			'sr_are': results['SR-ARE'],
			'sr_atad5': results["SR-ATAD5"],
			'sr_hse': results['SR-HSE'],
			'sr_mmp': results["SR-MMP"],
			'sr_p53': results["SR-p53"],
			'ld50': '%.3f'%results["LD50"],
			'similarity': '%.3f'%(sim*100.0),
			'pic': pic_url
		}
		return render(request, 'result.html', context=context)

