from genome_editing.flask_server.app import db


class SgrnaDesign(db.Model):
    __tablename__ = 'hg19_whole_genome_sgrna_design'
    index = db.Column(db.Integer, primary_key=True, unique=True)
    gene_symbol = db.Column(db.String(20))
    refseq_id = db.Column(db.String(50))
    exon_id = db.Column(db.Integer)
    chrom = db.Column(db.String(10))
    start = db.Column(db.Integer)
    end = db.Column(db.Integer)
    raw_sequence = db.Column(db.String(50))
    pam_type = db.Column(db.String(20))
    cutting_site_type = db.Column(db.String(20))
    cutting_site = db.Column(db.Numeric)
    sgrna_seq = db.Column(db.String(50))
    sgrna_full_seq = db.Column(db.String(50))
    sgrna_id = db.Column(db.Integer)
