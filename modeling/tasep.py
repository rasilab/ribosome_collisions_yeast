#!/usr/bin/env python3
"""Definition of Kinetic Model of mRNA Translation

Depends on rasilab customized versions of:
    - PySB (Model specification)
    - BioNetGen (Model parsing)
    - NFsim (Model simulation)

Install above from https://github.com/rasilab 

Arvind Rasi Subramaniam <rasi@fredhutch.org>
11 Sep 2018
"""

import pysb as sb

class Tasep(sb.Model):
    """Class to model motion of ribosomes on mRNAs.
    Derived from the generic pysb.Model

    Functions are defined for each type of reaction so that
    they can be reused with different parameter values.

    The 'tag' keyword for reaction definitions was added to PySB by Rasi to
    allow specific reactions to be monitored during the course of the
    simulation.

    Args:
        kwargs - dict of any/all parameters in the model below
    """

    def _set_params_(self):
        """Sets all parameters of the model

        All rates are specified in units of s-1 unless noted otherwise.
        All distances are specified in units of codons unless noted otherwise.
        """
        # molecule properties
        # length of ribosome
        sb.Parameter('l_ribo', 10)
        # length of mRNA
        sb.Parameter('l_mrna', 650)

        # canonical translation params
        # initiation rate
        sb.Parameter('k_init', 1)
        # elongation rate
        sb.Parameter('k_elong', 10)
        # termination rate
        sb.Parameter('k_term', 1)

        # canonical mrna metabolism params
        # transcription rate
        # rate of transcription is not important since it just controls the
        # number of mRNAs at any given time
        sb.Parameter('k_transcription', 0.001)
        # number of As in the polyA tail
        sb.Parameter('l_polya', 60)
        # rate at which each A in the polyA tail is removed
        sb.Parameter('k_deadenylation', 0.03)
        # rate of decapping once the polyA tail is completely removed
        sb.Parameter('k_decapping', 0.01)
        # rate at which each nt in the mRNA is removed 5' to 3'
        sb.Parameter('k_exo_53', 1)
        # rate at which each nt in the mRNA is removed 3' to 5'
        sb.Parameter('k_exo_35', 0)

        # stalling params
        # rate of elongation at the first stall
        # note that more stalls can be added by redefining the class
        sb.Parameter('k_elong_stall_1', 0.1)
        # location of first stall
        sb.Parameter('x_stall_1', 400)
        # number of stalls
        sb.Parameter('n_stall', 1)

        # premature termination params
        # rates of premature termination at intact A-sites
        # rate at which ribosomes that are not hit on either 5' or 3' abort
        sb.Parameter('k_preterm_no_hit_intact', 0.01)
        # rate at which ribosomes that are hit 5' abort
        sb.Parameter('k_preterm_5_hit_intact', 0.01)
        # rate at which ribosomes that are hit 3' abort
        sb.Parameter('k_preterm_3_hit_intact', 0.01)
        # rate at which ribosomes that are hit both 5' and 3' abort
        sb.Parameter('k_preterm_both_hit_intact', 0.01)
        # rates of premature termination at endocleaved A-sites
        # rate at which ribosomes that are not hit on either 5' or 3' abort
        sb.Parameter('k_preterm_no_hit_endocleaved', 1)
        # rate at which ribosomes that are hit 5' abort
        sb.Parameter('k_preterm_5_hit_endocleaved', 1)
        # rate at which ribosomes that are hit 3' abort
        sb.Parameter('k_preterm_3_hit_endocleaved', 1)
        # rate at which ribosomes that are hit both 5' and 3' abort
        sb.Parameter('k_preterm_both_hit_endocleaved', 1)

        # co-translational cleavage params
        # endocleavage occurs these many codons behind a-site
        sb.Parameter('l_cleave', 5)
        # rate at which ribosomes that are not hit on either 5' or 3' cause
        # cleavage
        sb.Parameter('k_cleave_no_hit', 0.0001)
        # rate at which ribosomes that are hit 5' cause cleavage
        sb.Parameter('k_cleave_5_hit', 0.0001)
        # rate at which ribosomes that are hit 3' cause cleavage
        sb.Parameter('k_cleave_3_hit', 0.0001)
        # rate at which ribosomes that are hit both 5' and 3' cause cleavage
        sb.Parameter('k_cleave_both_hit', 0.0001)

        # initial conditions
        # there is a single DNA locus
        sb.Parameter('n_dna_0', 1)
        # there are no mRNAs to begin with
        sb.Parameter('n_mrna_0', 0)
        # these many ribosomes
        # not critical since we set the initiation rate by hand above but keep
        # it high enough that there are enough ribosomes to translate all mRNAs
        # in the simulation
        sb.Parameter('n_ribosome_0', 10000)
        # no full proteins at the beginning
        sb.Parameter('n_protein_0', 0)
        # no aborted proteins at the beginning
        sb.Parameter('n_abortedprotein_0', 0)

    def _define_molecules_(self):
        """Define the molecules and their internal states
        """

        # rename first for nicer code
        self.ribosome_footprint_size = int(self.parameters['l_ribo'].value)
        self.cleave_distance = int(self.parameters['l_cleave'].value)
        self.mrna_length = int(self.parameters['l_mrna'].value)
        self.polya_length = int(self.parameters['l_polya'].value)
        elong_rates = [
            (int(self.parameters['x_stall_' + str(n + 1)].value),
             self.parameters['k_elong_stall_' + str(n + 1)])
            for n in range(int(self.parameters['n_stall'].value))]

        self.stall_elongation_rates = {k: v for (k, v) in elong_rates}

        # gene
        # the gene has no internal states, just used for transcription
        sb.Monomer('dna')

        # mRNA
        # specify sites and their states
        # start_region indicates whether mrna is initiable or not
        # the remaining mrna_length sites are codons for indicating
        # 1. 'p' is the asite and indicates ribosome occupancy
        # 2. 'r' is the RNA backbone and indicates cut or not cut
        #     cut can be caused by endo or exonucleolysis
        mrna_sites = ['cap', 'start_region']
        for pos in range(self.mrna_length):
            mrna_sites.append(f'p{pos}')
            mrna_sites.append(f'r{pos}')
        for pos in range(self.polya_length):
            mrna_sites.append(f'pA{pos}')

        mrna_site_states = {
            # start can be available to initiate or blocked
            'start_region': ['free', 'blocked'],
            # cap is necessary for initiation
            'cap': ['yes', 'no']
        }

        # this is to indicate if the mrna position has been inactivated by
        # endo/exo nucleolysis.
        for pos in range(self.mrna_length):
            mrna_site_states.update(
                {f'r{pos}': ['intact', 'endocleaved', 'exocleaved']}
            )
        # this is to indicate if the polyA position has been inactivated by
        # deadenylation.
        for pos in range(self.polya_length):
            mrna_site_states.update({f'pA{pos}' : ['intact', 'exocleaved']})

        # create the mRNA molecule
        sb.Monomer('mrna', sites=mrna_sites, site_states=mrna_site_states)

        # ribosome
        # specify sites and their states
        # the asite can be on mRNA or not
        # hit5 and hit3 indicate whether hit by another ribosome
        # on 5' or 3' resp.
        ribosome_sites = ['asite', 'hit3', 'hit5']
        # I am making hit5 an integer to remember how many times a ribosome has
        # been hit
        sb.Monomer('ribosome', sites=ribosome_sites)

        # the proteins have no internal states
        # full protein
        sb.Monomer('protein')
        # aborted protein
        sb.Monomer('abortedprotein')

    def _define_observables_(self):
        """Define the observables that are tracked duringn the simulation.

        Not critical for Rasi's version since we track all reactions
        """

        sb.Observable('n_protein', protein())
        sb.Observable('n_abortedprotein', abortedprotein())
        sb.Observable('n_ribosome_free', ribosome(asite=None))
        sb.Observable('n_ribosome', ribosome())
        sb.Observable('n_mrna', mrna())
        # track how many mrnas are inactive
        # (cannot produce proteins anymore once last codon is degraded)
        dead_mrna_args = {'r' + str(self.mrna_length - 1): 'exocleaved'}
        sb.Observable('n_mrna_dead', mrna(**dead_mrna_args))
        sb.Observable('n_ribosome_collision', ribosome(hit3=1) %
                      ribosome(hit5=1))

    def _define_initial_conditions_(self):
        """Define the initial number and states of molecules
        """
        sb.Initial(dna(), n_dna_0)

        initial_mrna = {
            'start_region': 'free',
            'cap': 'yes',
        }

        for pos in range(self.mrna_length):
            initial_mrna['r' + str(pos)] = 'intact'
            initial_mrna['p' + str(pos)] = None
        for pos in range(self.polya_length):
            initial_mrna['pA' + str(pos)] = 'intact'

        if self.parameters['n_mrna_0'].value > 0:
            sb.Initial(mrna(**initial_mrna), n_mrna_0)

        if self.parameters['n_protein_0'].value > 0:
            sb.Initial(protein(), n_protein_0)

        if self.parameters['n_abortedprotein_0'].value > 0:
            sb.Initial(abortedprotein(), n_abortedprotein_0)

        # ribosomes are intially free and not hit by other ribosomes
        initial_ribosome = {'asite': None, 'hit3': None, 'hit5': None}
        sb.Initial(ribosome(**initial_ribosome), n_ribosome_0)

    def transcribe(self, dna, mrna, k,
                   l_mrna=0, l_polya=0, tag=False):
        """Produce an mRNA from DNA
        """
        mrna_product_args = {'cap': 'yes',
                             'start_region': 'free'}

        # new mrnas do not have any ribosomes;
        # no sites either in coding or polyatail are cleaved
        for pos in range(l_mrna):
            mrna_product_args[f'p{pos}'] = None
            mrna_product_args[f'r{pos}'] = 'intact'
        for pos in range(l_polya):
            mrna_product_args['pA' + str(pos)] = 'intact'

        sb.Rule('transcription',
                dna() >> dna() + mrna(**mrna_product_args),
                k,
                tag=tag
                )

    def initiate(self, ribosome, mrna, pos, k, tag=False):
        """Initiate a ribosome at codon 'pos' on mRNA with rate k
        """
        # ribosomes can initiate only if the initation codon
        # has not been inactivated by cleavage
        mrna_reactant_args = {'p' + str(pos): None,
                              'r' + str(pos): 'intact',
                              'start_region': 'free',
                              'cap': 'yes'}
        mrna_product_args = {'p' + str(pos): 1,
                             'r' + str(pos): 'intact',
                             'start_region': 'blocked',
                             'cap': 'yes'}
        ribosome_reactant_args = {'asite': None}
        ribosome_product_args = {'asite': 1}
        sb.Rule(
            'initiation',
            ribosome(**ribosome_reactant_args) +
            mrna(**mrna_reactant_args) >>
            ribosome(**ribosome_product_args) %
            mrna(**mrna_product_args),
            k, total_rate=True, tag=tag)

    def elongate(self, ribosome, mrna, pos, k,
                 mrna_length, ribosome_footprint_size, tag=False):
        """Move ribosome from codon 'pos' on mRNA to 'pos' + 1
        """
        # ribosomes can elongate only if mrna is not inactivated by cleavage
        mrna_reactant_args = {
            f'p{pos}' : 1, f'p{pos + 1}': None, f'r{pos}': 'intact',
        }
        mrna_product_args = {
            f'p{pos}' : None, f'p{pos + 1}': 1, f'r{pos}': 'intact',
        }
        ribosome_reactant_args = {'asite': 1, 'hit5': None, 'hit3': None}
        ribosome_product_args = {'asite': 1, 'hit5': None, 'hit3': None}

        # ribosomes can elongate only if there is no ribosome immediately in
        # front this extra condition is necessary only if the ribosome is
        # larger than footprint
        if (pos < mrna_length - ribosome_footprint_size and
                ribosome_footprint_size > 1):
            mrna_reactant_args['p' + str(pos + ribosome_footprint_size)] = None
            mrna_product_args['p' + str(pos + ribosome_footprint_size)] = None

        # if ribosome moves beyond a footprint from ATG, mrna is free to
        # initiate
        if pos == ribosome_footprint_size - 1:
            mrna_reactant_args['start_region'] = 'blocked'
            mrna_product_args['start_region'] = 'free'
        sb.Rule(
            'elongation_' + str(pos),
            ribosome(**ribosome_reactant_args) %
            mrna(**mrna_reactant_args) >>
            ribosome(**ribosome_product_args) %
            mrna(**mrna_product_args),
            k, tag=tag)
        # elongation out of a collision, need to consider the ribosome behind
        if pos > ribosome_footprint_size - 1:
            mrna_reactant_args['p' +
                               str(pos - ribosome_footprint_size)] = 2
            mrna_product_args['p' +
                              str(pos - ribosome_footprint_size)] = 2
            ribosome_reactant_args = {
                'asite': 1,
                'hit5': 3,
                'hit3': None
            }
            ribosome_product_args = {
                'asite': 1,
                'hit5': None,
                'hit3': None
            }
            ribosome_5_reactant_args = {
                'asite': 2,
                'hit3': 3,
            }
            ribosome_5_product_args = {
                'asite': 2,
                'hit3': None,
            }
            sb.Rule(
                'elongation_with_hit5_' + str(pos),
                ribosome(**ribosome_5_reactant_args) %
                ribosome(**ribosome_reactant_args) %
                mrna(**mrna_reactant_args) >>
                ribosome(**ribosome_5_product_args) %
                ribosome(**ribosome_product_args) %
                mrna(**mrna_product_args),
                k, tag=tag)

    def terminate(self, ribosome, mrna, pos, k,
                  mrna_length, ribosome_footprint_size,
                  cleavestate,
                  tag=False):
        """Terminate ribosome by leaving mRNA at position pos.
        Ribosome can be not hit or hit from 5'.
        """
        mrna_reactant_args = {
            'p' + str(pos): 1,
        }
        mrna_product_args = {
            'p' + str(pos): None,
        }
        ribosome_reactant_args = {'asite': 1, 'hit3': None, 'hit5': None}
        ribosome_product_args = {'asite': None, 'hit3': None, 'hit5': None}
        # mrna becomes free to initiate if ribosomes terminated
        # within footprint distance of ATG
        if pos <= ribosome_footprint_size:
            mrna_reactant_args['start_region'] = 'blocked'
            mrna_product_args['start_region'] = 'free'
        # termination of ribosomes that are hit neither front or back
        sb.Rule(
            'term_no_hit_{}'.format(pos),
            ribosome(**ribosome_reactant_args) %
            mrna(**mrna_reactant_args) >>
            ribosome(**ribosome_product_args) +
            mrna(**mrna_product_args) +
            protein(),
            k,
            tag=tag)
        # termination out of a collision, need to consider the ribosome behind
        # but only if not separate by cleavage
        if pos > ribosome_footprint_size:
            mrna_reactant_args = {
                'p' + str(pos - ribosome_footprint_size): 1,
                'p' + str(pos): 3,
                'r' + str(pos): cleavestate,
            }
            mrna_product_args = {
                'p' + str(pos - ribosome_footprint_size): 1,
                'p' + str(pos): None,
                'r' + str(pos): cleavestate,
            }
            ribosome_reactant_args = {'asite': 3, 'hit3': None, 'hit5': 2}
            ribosome_product_args = {'asite': None, 'hit3': None, 'hit5': None}
            ribosome_5_reactant_args = {'asite': 1, 'hit3': 2}
            ribosome_5_product_args = {'asite': 1, 'hit3': None}
            sb.Rule(
                'term_5_hit_{}_{}'.format(cleavestate, pos),
                ribosome(**ribosome_5_reactant_args) %
                mrna(**mrna_reactant_args) %
                ribosome(**ribosome_reactant_args) >>
                ribosome(**ribosome_5_product_args) %
                mrna(**mrna_product_args) +
                ribosome(**ribosome_product_args) +
                protein(),
                k,
                tag=tag)

    def ribosomes_collide(self, ribosome, mrna, pos, k,
                          mrna_length, ribosome_footprint_size,
                          tag=False):
        """Ribosome at codon pos collides with ribosome at codon
        pos + ribosome_footprint_size on mRNA
        """
        if pos > mrna_length - ribosome_footprint_size:
            raise ValueError('There cannot be ribosome after stop!')
        # ribosomes can collide only if the back ribosome has an intact A-site
        mrna_reactant_args = {
            f'p{pos}': 1,
            f'r{pos}': 'intact',
            f'p{pos + ribosome_footprint_size}': 2,
        }
        mrna_product_args = {
            f'p{pos}': 1,
            f'r{pos}': 'intact',
            f'p{pos + ribosome_footprint_size}': 2,
        }
        ribosome_reactant_args = {'asite': 1, 'hit3': None}
        ribosome_product_args = {'asite': 1, 'hit3': 3}
        ribosome_3_reactant_args = {'asite': 2, 'hit5': None}
        ribosome_3_product_args = {'asite': 2, 'hit5': 3}
        sb.Rule(
            'collision_' + str(pos),
            ribosome(**ribosome_reactant_args) %
            ribosome(**ribosome_3_reactant_args) %
            mrna(**mrna_reactant_args) >>
            ribosome(**ribosome_product_args) %
            ribosome(**ribosome_3_product_args) %
            mrna(**mrna_product_args),
            k, tag=tag)

    def preterm_no_hit(self, ribosome, mrna, pos, k,
                       mrna_length, ribosome_footprint_size,
                       cleavestate, tag=False):
        """Premature termination of ribosomes that are hit neither
        from front or back.
        cleavestate indicate whether mrna pos is intact or cleaved.
        """
        mrna_reactant_args = {
            'p' + str(pos): 1,
            'r' + str(pos): cleavestate
        }
        mrna_product_args = {
            'p' + str(pos): None,
            'r' + str(pos): cleavestate
        }
        ribosome_reactant_args = {'asite': 1, 'hit3': None, 'hit5': None}
        ribosome_product_args = {'asite': None, 'hit3': None, 'hit5': None}
        # mrna becomes free to initiate if ribosomes terminated
        # within footprint distance of ATG
        if pos <= ribosome_footprint_size:
            mrna_reactant_args['start_region'] = 'blocked'
            mrna_product_args['start_region'] = 'free'
        # termination of ribosomes that are hit neither front or back
        sb.Rule(
            'preterm_no_hit_{}_{}'.format(cleavestate, pos),
            ribosome(**ribosome_reactant_args) %
            mrna(**mrna_reactant_args) >>
            ribosome(**ribosome_product_args) +
            mrna(**mrna_product_args) +
            abortedprotein(),
            k, tag=tag)

    def preterm_5_hit(self, ribosome, mrna, pos, k,
                      mrna_length, ribosome_footprint_size,
                      cleavestate, tag=False):
        """Premature termination of ribosomes that are hit only from 5'.
        cleavestate indicate whether mrna pos is intact or cleaved.
        """
        # note that mRNA cannot become free from this type of termination
        if pos < ribosome_footprint_size:
            raise ValueError("There cannot be ribosome behind this one!")
        mrna_reactant_args = {
            'p' + str(pos - ribosome_footprint_size): 1,
            'p' + str(pos): 3,
            'r' + str(pos): cleavestate,
        }
        mrna_product_args = {
            'p' + str(pos - ribosome_footprint_size): 1,
            'p' + str(pos): None,
            'r' + str(pos): cleavestate,
        }
        ribosome_reactant_args = {'asite': 3, 'hit3': None, 'hit5': 2}
        ribosome_product_args = {'asite': None, 'hit3': None, 'hit5': None}
        ribosome_5_reactant_args = {'asite': 1, 'hit3': 2}
        ribosome_5_product_args = {'asite': 1, 'hit3': None}
        sb.Rule(
            'preterm_5_hit_{}_{}'.format(cleavestate, pos),
            ribosome(**ribosome_5_reactant_args) %
            mrna(**mrna_reactant_args) %
            ribosome(**ribosome_reactant_args) >>
            ribosome(**ribosome_5_product_args) %
            mrna(**mrna_product_args) +
            ribosome(**ribosome_product_args) +
            abortedprotein(),
            k, tag=tag)

    def preterm_3_hit(self, ribosome, mrna, pos, k,
                      mrna_length, ribosome_footprint_size,
                      cleavestate, tag=False):
        """Premature termination of ribosomes that are hit only from 3'.
        cleavestate indicate whether mrna pos is intact or cleaved.
        """
        if pos >= mrna_length - ribosome_footprint_size:
            raise ValueError("There can be ribosome after stop!")
        mrna_reactant_args = {
            'p' + str(pos): 1,
            'r' + str(pos): cleavestate,
            'p' + str(pos + ribosome_footprint_size): 3,
        }

        mrna_product_args = {
            'p' + str(pos): None,
            'r' + str(pos): cleavestate,
            'p' + str(pos + ribosome_footprint_size): 3,
        }
        ribosome_reactant_args = {'asite': 1, 'hit3': 2, 'hit5': None}
        ribosome_product_args = {'asite': None, 'hit3': None, 'hit5': None}
        ribosome_3_reactant_args = {'asite': 3, 'hit5': 2}
        ribosome_3_product_args = {'asite': 3, 'hit5': None}
        # mrna becomes free to initiate if ribosomes terminated
        # within footprint distance of ATG
        if pos <= ribosome_footprint_size:
            mrna_reactant_args['start_region'] = 'blocked'
            mrna_product_args['start_region'] = 'free'
        sb.Rule(
            'preterm_3_hit_{}_{}'.format(cleavestate, pos),
            ribosome(**ribosome_reactant_args) %
            mrna(**mrna_reactant_args) %
            ribosome(**ribosome_3_reactant_args) >>
            ribosome(**ribosome_product_args) +
            mrna(**mrna_product_args) %
            ribosome(**ribosome_3_product_args) +
            abortedprotein(),
            k, tag=tag)

    def preterm_both_hit(self, ribosome, mrna, pos, k,
                         mrna_length, ribosome_footprint_size,
                         cleavestate, tag=False):
        """Premature termination of ribosomes that are hit both 5' and 3'.
        cleavestate indicate whether mrna pos is intact or cleaved.
        """
        if pos >= mrna_length - ribosome_footprint_size:
            raise ValueError("There can be ribosome after stop!")

        if pos < ribosome_footprint_size:
            raise ValueError("There cannot be ribosome behind this one!")

        mrna_reactant_args = {
            'p' + str(pos - ribosome_footprint_size): 1,
            'p' + str(pos): 3,
            'r' + str(pos): cleavestate,
            'p' + str(pos + ribosome_footprint_size): 5,
        }

        mrna_product_args = {
            'p' + str(pos - ribosome_footprint_size): 1,
            'p' + str(pos): None,
            'r' + str(pos): cleavestate,
            'p' + str(pos + ribosome_footprint_size): 5,
        }

        ribosome_5_reactant_args = {'asite': 1, 'hit3': 2}
        ribosome_5_product_args = {'asite': 1, 'hit3': None}
        ribosome_reactant_args = {'asite': 3, 'hit3': 4, 'hit5': 2}
        ribosome_product_args = {'asite': None, 'hit3': None, 'hit5': None}
        ribosome_3_reactant_args = {'asite': 5, 'hit5': 4}
        ribosome_3_product_args = {'asite': 5, 'hit5': None}

        sb.Rule(
            'preterm_both_hit_{}_{}'.format(cleavestate, pos),
            ribosome(**ribosome_5_reactant_args) %
            ribosome(**ribosome_reactant_args) %
            mrna(**mrna_reactant_args) %
            ribosome(**ribosome_3_reactant_args) >>
            ribosome(**ribosome_product_args) +
            ribosome(**ribosome_5_product_args) %
            mrna(**mrna_product_args) %
            ribosome(**ribosome_3_product_args) +
            abortedprotein(),
            k, tag=tag)

    def endocleave_no_hit(self, ribosome, mrna, pos, k,
                          mrna_length, ribosome_footprint_size, cleave_distance,
                          tag=False):
        """Cleave mRNA behind the ribosome which is not hit.
        Note that we assume that only one endonucleolytic cleavage
        occurs per mRNA. So we model this as only capped mRNAs
        can be endonucleolytically cleaved.
        The mRNA also becomes decapped right away so that
        it cannot be initiated.
        """
        mrna_reactant_args = {
            'p' + str(pos): 1,
            'r' + str(pos - ribosome_footprint_size): 'intact',
            'cap': 'yes',
        }
        mrna_product_args = {
            'p' + str(pos): 1,
            'r' + str(pos - ribosome_footprint_size): 'endocleaved',
            'cap': 'no'
        }
        ribosome_reactant_args = {'asite': 1, 'hit3': None, 'hit5': None}
        ribosome_product_args = {'asite': 1, 'hit3': None, 'hit5': None}
        sb.Rule(
            'endocleave_no_hit_' + str(pos),
            ribosome(**ribosome_reactant_args) %
            mrna(**mrna_reactant_args) >>
            ribosome(**ribosome_product_args) %
            mrna(**mrna_product_args),
            k, tag=tag)

    def endocleave_5_hit(self, ribosome, mrna, pos, k,
                         mrna_length, ribosome_footprint_size, cleave_distance,
                         tag=False):
        """Cleave mRNA behind the ribosome which is 5' hit.
        Bond between the 5' ribosome and this ribosome is broken.
        Note that we assume that only one endonucleolytic cleavage
        occurs per mRNA. So we model this as only capped mRNAs
        can be endonucleolytically cleaved.
        The mRNA also becomes decapped right away so that
        it cannot be initiated.
        """
        mrna_reactant_args = {
            'p' + str(pos - ribosome_footprint_size): 1,
            'r' + str(pos - ribosome_footprint_size): 'intact',
            'p' + str(pos): 3,
            'cap': 'yes'
        }
        mrna_product_args = {
            'p' + str(pos - ribosome_footprint_size): 1,
            'r' + str(pos - ribosome_footprint_size): 'endocleaved',
            'p' + str(pos): 3,
            'cap': 'no'
        }
        ribosome_5_reactant_args = {'asite': 1, 'hit3': 2}
        ribosome_5_product_args = {'asite': 1, 'hit3': None}
        ribosome_reactant_args = {'asite': 3, 'hit5': 2, 'hit3': None}
        ribosome_product_args = {'asite': 3, 'hit5': None, 'hit3': None}
        sb.Rule(
            'endocleave_5_hit_{}'.format(pos),
            ribosome(**ribosome_5_reactant_args) %
            mrna(**mrna_reactant_args) %
            ribosome(**ribosome_reactant_args) >>
            ribosome(**ribosome_5_product_args) %
            mrna(**mrna_product_args) %
            ribosome(**ribosome_product_args),
            k, tag=tag)

    def endocleave_3_hit(self, ribosome, mrna, pos, k,
                         mrna_length, ribosome_footprint_size, cleave_distance,
                         tag=False):
        """Cleave mRNA behind the ribosome which is 3' hit.
        Note that we assume that only one endonucleolytic cleavage
        occurs per mRNA. So we model this as only capped mRNAs
        can be endonucleolytically cleaved.
        The mRNA also becomes decapped right away so that
        it cannot be initiated.
        """
        mrna_reactant_args = {
            'p' + str(pos): 1,
            'r' + str(pos - ribosome_footprint_size): 'intact',
            'p' + str(pos + ribosome_footprint_size): 3,
            'cap': 'yes'
        }
        mrna_product_args = {
            'p' + str(pos): 1,
            'r' + str(pos - ribosome_footprint_size): 'endocleaved',
            'p' + str(pos + ribosome_footprint_size): 3,
            'cap': 'no'
        }
        ribosome_reactant_args = {'asite': 1, 'hit5': None, 'hit3': 2}
        ribosome_product_args = {'asite': 1, 'hit5': None, 'hit3': 2}
        ribosome_3_reactant_args = {'asite': 3, 'hit5': 2}
        ribosome_3_product_args = {'asite': 3, 'hit5': 2}
        sb.Rule(
            'endocleave_3_hit_{}'.format(pos),
            ribosome(**ribosome_reactant_args) %
            mrna(**mrna_reactant_args) %
            ribosome(**ribosome_3_reactant_args) >>
            ribosome(**ribosome_product_args) %
            mrna(**mrna_product_args) %
            ribosome(**ribosome_3_product_args),
            k, tag=tag)

    def endocleave_both_hit(self, ribosome, mrna, pos, k,
                            mrna_length, ribosome_footprint_size, cleave_distance,
                            tag=False):
        """Cleave mRNA behind the ribosome which is hit both 5' and 3'.
        Bond between the 5' ribosome and this ribosome is broken.
        Note that we assume that only one endonucleolytic cleavage
        occurs per mRNA. So we model this as only capped mRNAs
        can be endonucleolytically cleaved.
        The mRNA also becomes decapped right away so that
        it cannot be initiated.
        """
        mrna_reactant_args = {
            'p' + str(pos - ribosome_footprint_size): 1,
            'r' + str(pos - ribosome_footprint_size): 'intact',
            'p' + str(pos): 3,
            'p' + str(pos + ribosome_footprint_size): 5,
            'cap': 'yes'
        }
        mrna_product_args = {
            'p' + str(pos - ribosome_footprint_size): 1,
            'r' + str(pos - ribosome_footprint_size): 'endocleaved',
            'p' + str(pos): 3,
            'p' + str(pos + ribosome_footprint_size): 5,
            'cap': 'no'
        }
        ribosome_5_reactant_args = {'asite': 1, 'hit3': 2}
        ribosome_5_product_args = {'asite': 1, 'hit3': None}
        ribosome_reactant_args = {'asite': 3, 'hit5': 2, 'hit3': 4}
        ribosome_product_args = {'asite': 3, 'hit5': None, 'hit3': 4}
        ribosome_3_reactant_args = {'asite': 5, 'hit5': 4}
        ribosome_3_product_args = {'asite': 5, 'hit5': 4}
        sb.Rule(
            'endocleave_both_hit_{}'.format(pos),
            ribosome(**ribosome_5_reactant_args) %
            ribosome(**ribosome_reactant_args) %
            mrna(**mrna_reactant_args) %
            ribosome(**ribosome_3_reactant_args) >>
            ribosome(**ribosome_5_product_args) %
            ribosome(**ribosome_product_args) %
            mrna(**mrna_product_args) %
            ribosome(**ribosome_3_product_args),
            k, tag=tag)

    def deadenylate(self, mrna, pos, k, l_polya, tag=False):
        """Inactivate polyA tail at position pos
        """
        # convert the position from intact to exocleaved
        mrna_reactant_args = {'pA' + str(pos): 'intact'}
        mrna_product_args = {'pA' + str(pos): 'exocleaved'}
        # deadenylation at non-terminal As can occur only
        # if previous A was deadenylated
        if pos < l_polya - 1:
            mrna_reactant_args['pA' + str(pos + 1)] = 'exocleaved'
            mrna_product_args['pA' + str(pos + 1)] = 'exocleaved'
        sb.Rule(
            'deadenylation_' + str(pos),
            mrna(**mrna_reactant_args) >> mrna(**mrna_product_args),
            k, tag=tag
        )

    def decap(self, mrna, k, tag=False):
        """Decap deadenylated mRNA
        """
        # decapping can occur only if the polyA at pos was deadenylated
        pos = 0
        mrna_reactant_args = {'pA' + str(pos): 'exocleaved'}
        mrna_product_args = {'pA' + str(pos): 'exocleaved'}
        mrna_reactant_args['cap'] = 'yes'
        mrna_product_args['cap'] = 'no'
        sb.Rule(
            'decapping',
            mrna(**mrna_reactant_args) >> mrna(**mrna_product_args),
            k, tag=tag
        )

    def exo_53(self, mrna, pos, k, tag=False):
        """Inactivate decapped or endonucleolytically mRNA from 5' to 3'
        """
        # convert the position from intact to exocleaved only if not occupied by
        # ribosome
        mrna_reactant_args = {f'p{pos}': None, f'r{pos}': 'intact'}
        mrna_product_args = {f'p{pos}': None, f'r{pos}': 'exocleaved'}
        # 5' to 3' exonucleolysis can occur only
        # if previous nt was also cleaved endo or exonucleolytically.
        # The exception is the first nt of the mRNA which can be
        # exonucleosed only if the cap is absent.
        if pos == 0:
            mrna_reactant_args['cap'] = 'no'
            mrna_product_args['cap'] = 'no'
            sb.Rule(
                'exonucleolysis_5end_' + str(pos),
                mrna(**mrna_reactant_args) >> mrna(**mrna_product_args),
                k, tag=tag
            )
        else:
            # exo
            mrna_reactant_args[f'r{pos - 1}'] = 'exocleaved'
            mrna_product_args[f'r{pos - 1}'] = 'exocleaved'
            sb.Rule(
                'exonucleolysis_5end_' + str(pos),
                mrna(**mrna_reactant_args) >> mrna(**mrna_product_args),
                k, tag=tag
            )
            # endo
            mrna_reactant_args[f'r{pos - 1}'] = 'endocleaved'
            mrna_product_args[f'r{pos - 1}'] = 'endocleaved'
            sb.Rule(
                'exonucleolysis_from_endo_5end_' + str(pos),
                mrna(**mrna_reactant_args) >> mrna(**mrna_product_args),
                k, tag=tag
            )

    def exo_35(self, mrna, pos, k, tag=False):
        """Inactivate deadenylated/endonucleolytically cleaved mRNA from 3' to 5'
        """
        # convert the position from intact to exocleaved only if not occupied by
        # ribosome
        mrna_reactant_args = {f'p{pos}': None, f'r{pos}': 'intact'}
        mrna_product_args = {f'p{pos}': None, f'r{pos}': 'exocleaved'}
        # 3' to 5' exonucleolysis can occur only
        # if 3' nt was also cleaved endo or exonucleolytically.
        # The exception is the last nt of the mRNA which can be
        # exonucleosed only if the poly A is absent.
        if pos == self.mrna_length:
            mrna_reactant_args['pA1'] = 'exocleaved'
            mrna_product_args['pA1'] = 'exocleaved'
            sb.Rule(
                'exonucleolysis_3end_' + str(pos),
                mrna(**mrna_reactant_args) >> mrna(**mrna_product_args),
                k, tag=tag
            )
        else:
            # exo
            mrna_reactant_args[f'r{pos + 1}'] = 'exocleaved'
            mrna_product_args[f'r{pos + 1}'] = 'exocleaved'
            sb.Rule(
                'exonucleolysis_3end_' + str(pos),
                mrna(**mrna_reactant_args) >> mrna(**mrna_product_args),
                k, tag=tag
            )
            # endo
            mrna_reactant_args[f'r{pos + 1}'] = 'endocleaved'
            mrna_product_args[f'r{pos + 1}'] = 'endocleaved'
            sb.Rule(
                'exonucleolysis_from_endo_3end_' + str(pos),
                mrna(**mrna_reactant_args) >> mrna(**mrna_product_args),
                k, tag=tag
            )

    def _redefine_params_(self, kwargs):
        """Redefine params using dictionary passed during model instantiation.
        """
        # reset parameters if passed as keyword during model instantiation
        for k, v in kwargs.items():
            if k not in ['x_stall', 'k_elong_stall']:
                try:
                    self.parameters[k].value = v
                except KeyError:
                    continue
            elif k == 'x_stall':
                if isinstance(v, str):
                    x_stall = [int(x) for x in v.split(',')]
                elif isinstance(v, float) or isinstance(v, int):
                    x_stall = [int(v)]
                else:
                    raise(ValueError("x_stall must be an integer or "
                                     "comma-separated string of integers"))
                # reset the first stall location
                self.parameters['x_stall_1'].value = x_stall[0]
                # create one param for each additional stall location
                for n in range(1, len(x_stall)):
                    sb.Parameter(f'x_stall_{n+1}', x_stall[n])
                self.parameters['n_stall'].value = len(x_stall)
            elif k == 'k_elong_stall':
                if isinstance(v, str):
                    k_elong_stall = [float(x) for x in v.split(',')]
                elif isinstance(v, float) or isinstance(v, int):
                    k_elong_stall = [v]
                else:
                    raise(ValueError("k_elong_stall must be an integer or "
                                     "comma-separated string of integers"))
                # reset the first stall duration 
                self.parameters['k_elong_stall_1'].value = k_elong_stall[0]
                # create one param for each additional stall location
                for n in range(1, len(k_elong_stall)):
                    sb.Parameter(f'k_elong_stall_{n + 1}', k_elong_stall[n])

    def __init__(self, **kwargs):
        """Instantiate an object of the class by calling above functions

        We use loops to define reactions that occur at multiple positions
        along mRNA.
        """

        # instantiate a pysb model object with default minimal structure
        super().__init__()

        # set parameters for the model to their default values
        self._set_params_()

        # redefine parameters if they were passed during instantiation.
        # this is often done to vary the parameters of the model
        self._redefine_params_(kwargs)

        # define molecules, their states and state values
        self._define_molecules_()

        # define observables -- states that are tracked during the simulation
        self._define_observables_()

        # define initial copy number of all molecular states
        self._define_initial_conditions_()

        ############################
        # Define Reactions / Rules #
        ############################

        # transcription
        self.transcribe(dna, mrna, k_transcription,
                        self.mrna_length, self.polya_length, tag=True)

        # translation initiation
        self.initiate(ribosome, mrna, 0, k_init, tag=True)

        # elongation
        for pos in range(self.mrna_length - 1):
            # normal
            if pos not in self.stall_elongation_rates:
                self.elongate(ribosome, mrna, pos, k_elong,
                              self.mrna_length, self.ribosome_footprint_size, tag=True)
            # stalled
            else:
                self.elongate(ribosome, mrna, pos, self.stall_elongation_rates[pos],
                              self.mrna_length, self.ribosome_footprint_size, tag=True)

        # normal termination from intact stop codons (not endocleaved)
        self.terminate(ribosome, mrna, self.mrna_length - 1, k_term,
                       self.mrna_length, self.ribosome_footprint_size,
                       "intact", tag=True)

        # deadenylation
        for pos in range(self.polya_length):
            self.deadenylate(mrna, pos, k_deadenylation,
                             self.polya_length, tag=True)

        # decapping
        self.decap(mrna, k_decapping, tag=True)

        # 5' to 3' exonucleolysis
        for pos in range(self.mrna_length):
            self.exo_53(mrna, pos, k_exo_53, tag=True)

        # 3' to 5' exonucleolysis
        for pos in range(self.mrna_length):
            self.exo_35(mrna, pos, k_exo_35, tag=True)

        # collision
        for pos in range(self.mrna_length - self.ribosome_footprint_size):
            # normal sites
            if pos not in self.stall_elongation_rates:
                self.ribosomes_collide(ribosome, mrna, pos, k_elong,
                                       self.mrna_length,
                                       self.ribosome_footprint_size,
                                       tag=True)
            # stall sites
            else:
                self.ribosomes_collide(ribosome, mrna, pos, self.stall_elongation_rates[pos],
                                       self.mrna_length,
                                       self.ribosome_footprint_size,
                                       tag=True)

        # abortive/premature termination of ribosomes that are not hit
        for pos in range(self.mrna_length):
            # intact mRNA at A-site
            self.preterm_no_hit(ribosome, mrna, pos, k_preterm_no_hit_intact,
                                self.mrna_length,
                                self.ribosome_footprint_size,
                                "intact", tag=True)
            # cleaved mRNA at A-site
            self.preterm_no_hit(ribosome, mrna, pos, k_preterm_no_hit_endocleaved,
                                self.mrna_length,
                                self.ribosome_footprint_size,
                                "endocleaved", tag=True)

        # cleavage of mrna at ribosomes that are not hit
        for pos in range(self.ribosome_footprint_size, self.mrna_length):
            self.endocleave_no_hit(ribosome, mrna, pos, k_cleave_no_hit,
                                   self.mrna_length,
                                   self.ribosome_footprint_size,
                                   self.cleave_distance, tag=True)

        # termination of ribosomes that are 5 hit
        for pos in range(self.ribosome_footprint_size, self.mrna_length):
            # intact mRNA at A-site
            self.preterm_5_hit(ribosome, mrna, pos, k_preterm_5_hit_intact,
                               self.mrna_length,
                               self.ribosome_footprint_size,
                               "intact", tag=True)
            # cleaved mRNA at A-site
            self.preterm_5_hit(ribosome, mrna, pos, k_preterm_5_hit_endocleaved,
                               self.mrna_length,
                               self.ribosome_footprint_size,
                               "endocleaved", tag=True)

        # cleavage of mrna at ribosomes that are 5 hit
        for pos in range(self.ribosome_footprint_size, self.mrna_length):
            self.endocleave_5_hit(ribosome, mrna, pos, k_cleave_5_hit,
                                  self.mrna_length,
                                  self.ribosome_footprint_size,
                                  self.cleave_distance, tag=True)

        # termination of ribosomes that are 3 hit
        for pos in range(self.mrna_length - self.ribosome_footprint_size):
            # intact mRNA at A-site
            self.preterm_3_hit(ribosome, mrna, pos, k_preterm_3_hit_intact,
                               self.mrna_length,
                               self.ribosome_footprint_size,
                               "intact", tag=True)
            # cleaved mRNA at A-site
            self.preterm_3_hit(ribosome, mrna, pos, k_preterm_3_hit_endocleaved,
                               self.mrna_length,
                               self.ribosome_footprint_size,
                               "endocleaved", tag=True)

        # cleavage of mrna at ribosomes that are 3 hit
        for pos in range(self.ribosome_footprint_size,
                         self.mrna_length - self.ribosome_footprint_size):
            self.endocleave_3_hit(ribosome, mrna, pos, k_cleave_3_hit,
                                  self.mrna_length,
                                  self.ribosome_footprint_size,
                                  self.cleave_distance, tag=True)

        # termination of ribosomes that are both hit
        for pos in range(self.ribosome_footprint_size,
                         self.mrna_length - self.ribosome_footprint_size):
            # intact mRNA at A-site
            self.preterm_both_hit(ribosome, mrna, pos, k_preterm_both_hit_intact,
                                  self.mrna_length,
                                  self.ribosome_footprint_size,
                                  "intact", tag=True)
            # cleaved mRNA at A-site
            self.preterm_both_hit(ribosome, mrna, pos, k_preterm_both_hit_endocleaved,
                                  self.mrna_length,
                                  self.ribosome_footprint_size,
                                  "endocleaved", tag=True)

        # cleavage of mrna at ribosomes that are both hit
        for pos in range(self.ribosome_footprint_size,
                         self.mrna_length - self.ribosome_footprint_size):
            self.endocleave_both_hit(ribosome, mrna, pos, k_cleave_both_hit,
                                     self.mrna_length,
                                     self.ribosome_footprint_size,
                                     self.cleave_distance, tag=True)

if __name__ == '__main__':
    """Run simulation with default values if the file is run as a script.

    Will create output files in the current working directory

    """

    import argparse
    import sys
    import subprocess as sp
    from pysb.export import export

    parser = argparse.ArgumentParser(description='Run translation simulation.')
    parser.add_argument('--bng2', type=str, metavar='BNG2.pl', 
                        help='Path to BioNetGen BNG2.pl script (default: BNG2.pl)',
                        default='BNG2.pl')
    parser.add_argument('--nfsim', type=str, metavar='NFsim',
                        help='Path to NFsim executable (default: NFsim)',
                        default='NFsim')
    args = parser.parse_args()

    equilibrium_time = 0  # seconds
    tstop = str(100000)  # seconds
    maxcputime = str(100 * 60)  # seconds
    osteps = str(10)  # number of samples
    seed = str(111)  # random number initial seed
    gml = str(1000000)  # max num of mol allowed in simulation
    utl = '3'  # max number of bonds to traverse during simulation
    network = '-connect'  # whether to infer reaction network connectivity

    # instantiate the model class with default parameter values, but change
    # mRNA length to 500 codons as an example of how to change parameters
    model = Tasep(**{'l_mrna': 500})

    outdir = '.'
    bnglfile = f'{outdir}/tasep.bngl'
    xmlfile = bnglfile.replace('.bngl', '.xml')
    gdatfile = bnglfile.replace('.bngl', '.gdat')
    rxnfile = bnglfile.replace('.bngl', '.rxns.tsv')
    paramsfile =  bnglfile.replace('.bngl', '.params.tsv')

    # write all model parameters to a separate file
    with open(paramsfile, 'w') as file:
        file.write('parameter\tvalue\n')
        for param in model.parameters:
            file.write(f'{param.name}\t{param.value}\n')
    # compress the params file
    sp.run(['gzip', '-f', paramsfile], cwd=outdir)

    # write BNGL file
    with open(bnglfile, 'w') as file:
        file.write(export(model, 'bngl'))

    # convert BNGL file to XML for NFSim input
    sp.run([args.bng2, '--xml', '--outdir', outdir, bnglfile])

    # simlate with NFSim
    # print NFSim command
    nfsim_command = [
        args.nfsim, '-xml', xmlfile, '-sim', tstop, '-oSteps', osteps,
        '-seed', seed, '-o', gdatfile, '-rxnlog', rxnfile,
        '-utl', utl,
        '-gml', gml, '-maxcputime', maxcputime,
        network
    ]
    print(" ".join(nfsim_command))

    # simlate with NFSim
    sp.run(nfsim_command)
